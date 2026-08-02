#include "ipc.h"

#include <chrono>
#include <condition_variable>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <thread>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/socket.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/un.h>
#include <unistd.h>
#include <cerrno>
#endif

namespace astap {
  namespace ipc {
    namespace {
      // A request is a file name and a handful of numbers, a reply a few lines
      // of text. The cap is here only so that a confused sender cannot make the
      // receiver allocate without bound.
      constexpr uint32_t kMaxMessage = 64u * 1024 * 1024;

      void encode_length(uint32_t n, unsigned char *out) {
        out[0] = static_cast<unsigned char>(n & 0xff);
        out[1] = static_cast<unsigned char>((n >> 8) & 0xff);
        out[2] = static_cast<unsigned char>((n >> 16) & 0xff);
        out[3] = static_cast<unsigned char>((n >> 24) & 0xff);
      }

      uint32_t decode_length(const unsigned char *in) {
        return static_cast<uint32_t>(in[0]) | (static_cast<uint32_t>(in[1]) << 8) |
               (static_cast<uint32_t>(in[2]) << 16) | (static_cast<uint32_t>(in[3]) << 24);
      }
    } // namespace

#ifdef _WIN32
    namespace {
      std::string last_error_text(DWORD e) {
        char *buf = nullptr;
        const DWORD n = FormatMessageA(
            FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM |
                FORMAT_MESSAGE_IGNORE_INSERTS,
            nullptr, e, 0, reinterpret_cast<char *>(&buf), 0, nullptr);
        std::string s = n && buf ? std::string(buf, n) : "error " + std::to_string(e);
        if (buf) LocalFree(buf);
        while (!s.empty() && (s.back() == '\n' || s.back() == '\r' || s.back() == '.')) s.pop_back();
        return s;
      }

      // Blocking I/O on an overlapped handle, with a deadline. The pipe is
      // opened overlapped purely so that a wedged peer costs a bounded wait
      // rather than a hung process; nothing here is asynchronous in spirit.
      bool overlapped_io(HANDLE h, void *buf, DWORD len, bool writing, int timeout_ms,
                         std::string *error) {
        HANDLE ev = CreateEventA(nullptr, TRUE, FALSE, nullptr);
        if (!ev) {
          if (error) *error = "CreateEvent: " + last_error_text(GetLastError());
          return false;
        }
        OVERLAPPED ov = {};
        ov.hEvent = ev;
        DWORD done = 0;
        BOOL ok = writing ? WriteFile(h, buf, len, &done, &ov) : ReadFile(h, buf, len, &done, &ov);
        if (!ok && GetLastError() == ERROR_IO_PENDING) {
          const DWORD w = WaitForSingleObject(ev, timeout_ms < 0 ? INFINITE
                                                                 : static_cast<DWORD>(timeout_ms));
          if (w != WAIT_OBJECT_0) {
            CancelIo(h);
            WaitForSingleObject(ev, INFINITE); // the cancel still completes it
            CloseHandle(ev);
            if (error) *error = w == WAIT_TIMEOUT ? "timed out" : "wait failed";
            return false;
          }
          ok = GetOverlappedResult(h, &ov, &done, FALSE);
        }
        const DWORD err = GetLastError();
        CloseHandle(ev);
        if (!ok) {
          if (error) *error = last_error_text(err);
          return false;
        }
        if (done != len) {
          if (error) *error = "short transfer";
          return false;
        }
        return true;
      }

      bool send_message(HANDLE h, const std::string &m, int timeout_ms, std::string *error) {
        unsigned char hdr[4];
        encode_length(static_cast<uint32_t>(m.size()), hdr);
        if (!overlapped_io(h, hdr, 4, true, timeout_ms, error)) return false;
        if (m.empty()) return true;
        return overlapped_io(h, const_cast<char *>(m.data()), static_cast<DWORD>(m.size()), true,
                             timeout_ms, error);
      }

      bool receive_message(HANDLE h, std::string &m, int timeout_ms, std::string *error) {
        unsigned char hdr[4];
        if (!overlapped_io(h, hdr, 4, false, timeout_ms, error)) return false;
        const uint32_t n = decode_length(hdr);
        if (n > kMaxMessage) {
          if (error) *error = "message too large";
          return false;
        }
        m.assign(n, '\0');
        if (!n) return true;
        return overlapped_io(h, &m[0], n, false, timeout_ms, error);
      }
    } // namespace

    std::string default_endpoint() { return "\\\\.\\pipe\\faster-astap-index"; }

    bool request(const std::string &endpoint, const std::string &message, std::string &reply,
                 int connect_timeout_ms, int reply_timeout_ms, std::string *error) {
      const auto deadline =
          std::chrono::steady_clock::now() + std::chrono::milliseconds(connect_timeout_ms);
      HANDLE h = INVALID_HANDLE_VALUE;
      for (;;) {
        h = CreateFileA(endpoint.c_str(), GENERIC_READ | GENERIC_WRITE, 0, nullptr, OPEN_EXISTING,
                        FILE_FLAG_OVERLAPPED, nullptr);
        if (h != INVALID_HANDLE_VALUE) break;
        const DWORD e = GetLastError();
        // BUSY means a server is there but every instance is serving someone
        // else; NOT_FOUND means there is no server yet, which is the case worth
        // waiting through while one starts up.
        if (e != ERROR_PIPE_BUSY && e != ERROR_FILE_NOT_FOUND) {
          if (error) *error = last_error_text(e);
          return false;
        }
        if (std::chrono::steady_clock::now() >= deadline) {
          if (error) *error = e == ERROR_PIPE_BUSY ? "server busy" : "no server listening";
          return false;
        }
        if (e == ERROR_PIPE_BUSY) WaitNamedPipeA(endpoint.c_str(), 50);
        else std::this_thread::sleep_for(std::chrono::milliseconds(20));
      }

      bool ok = send_message(h, message, reply_timeout_ms, error) &&
                receive_message(h, reply, reply_timeout_ms, error);
      CloseHandle(h);
      return ok;
    }

    bool endpoint_alive(const std::string &endpoint) {
      return WaitNamedPipeA(endpoint.c_str(), 0) != 0 || GetLastError() == ERROR_SEM_TIMEOUT;
    }

    struct Server::Impl {
      std::vector<HANDLE> pipes;
      std::vector<std::thread> workers;
      Handler handler;
      std::atomic<bool> stopping{false};
      std::mutex done_mutex;
      std::condition_variable done;

      void serve_one(HANDLE h) {
        std::string req;
        if (!receive_message(h, req, 300000, nullptr)) return;
        const std::string rep = handler(req);
        send_message(h, rep, 300000, nullptr);
        FlushFileBuffers(h);
      }
    };

    Server::Server() : impl_(new Impl) {}

    Server::~Server() { stop(); }

    bool Server::start(const std::string &endpoint, int instances, Handler handler,
                       std::string *error) {
      endpoint_ = endpoint;
      impl_->handler = std::move(handler);
      if (instances < 1) instances = 1;
      for (int i = 0; i < instances; i++) {
        DWORD open_mode = PIPE_ACCESS_DUPLEX | FILE_FLAG_OVERLAPPED;
        // Claiming the name on the first instance is what makes a second server
        // fail loudly here instead of quietly stealing half the clients.
        if (i == 0) open_mode |= FILE_FLAG_FIRST_PIPE_INSTANCE;
        HANDLE h = CreateNamedPipeA(endpoint.c_str(), open_mode,
                                    PIPE_TYPE_BYTE | PIPE_READMODE_BYTE | PIPE_WAIT,
                                    static_cast<DWORD>(instances), 64 * 1024, 64 * 1024, 0,
                                    nullptr);
        if (h == INVALID_HANDLE_VALUE) {
          const DWORD e = GetLastError();
          for (HANDLE p : impl_->pipes) CloseHandle(p);
          impl_->pipes.clear();
          if (error)
            *error = e == ERROR_ACCESS_DENIED
                         ? "another server already holds " + endpoint
                         : "CreateNamedPipe: " + last_error_text(e);
          return false;
        }
        impl_->pipes.push_back(h);
      }
      return true;
    }

    void Server::run() {
      for (HANDLE h : impl_->pipes) {
        impl_->workers.emplace_back([this, h] {
          while (!impl_->stopping) {
            OVERLAPPED ov = {};
            ov.hEvent = CreateEventA(nullptr, TRUE, FALSE, nullptr);
            BOOL ok = ConnectNamedPipe(h, &ov);
            DWORD e = GetLastError();
            if (!ok && e == ERROR_IO_PENDING) {
              WaitForSingleObject(ov.hEvent, INFINITE);
              DWORD n = 0;
              ok = GetOverlappedResult(h, &ov, &n, FALSE);
              e = GetLastError();
            }
            CloseHandle(ov.hEvent);
            // A client that connected in the window between creating the
            // instance and asking for a connection is already here, not an
            // error.
            if (ok || e == ERROR_PIPE_CONNECTED) {
              if (!impl_->stopping) impl_->serve_one(h);
              DisconnectNamedPipe(h);
            }
          }
        });
      }
      std::unique_lock<std::mutex> lk(impl_->done_mutex);
      impl_->done.wait(lk, [this] { return impl_->stopping.load(); });
      lk.unlock();

      // Every worker is parked in ConnectNamedPipe, and the only way back out is
      // for somebody to connect. One throwaway connection per instance does it.
      for (size_t i = 0; i < impl_->pipes.size(); i++) {
        HANDLE h = CreateFileA(endpoint_.c_str(), GENERIC_READ | GENERIC_WRITE, 0, nullptr,
                               OPEN_EXISTING, 0, nullptr);
        if (h != INVALID_HANDLE_VALUE) CloseHandle(h);
      }
      for (std::thread &t : impl_->workers)
        if (t.joinable()) t.join();
      impl_->workers.clear();
      for (HANDLE h : impl_->pipes) CloseHandle(h);
      impl_->pipes.clear();
    }

    void Server::stop() {
      if (impl_->stopping.exchange(true)) return;
      std::lock_guard<std::mutex> lk(impl_->done_mutex);
      impl_->done.notify_all();
    }

#else // --- POSIX -------------------------------------------------------------

    namespace {
      std::string errno_text() { return std::strerror(errno); }

      bool write_all(int fd, const char *p, size_t n) {
        while (n) {
          const ssize_t w = ::send(fd, p, n, 0);
          if (w <= 0) return false;
          p += w;
          n -= static_cast<size_t>(w);
        }
        return true;
      }

      bool read_all(int fd, char *p, size_t n) {
        while (n) {
          const ssize_t r = ::recv(fd, p, n, 0);
          if (r <= 0) return false;
          p += r;
          n -= static_cast<size_t>(r);
        }
        return true;
      }

      void set_timeout(int fd, int ms) {
        if (ms < 0) return;
        timeval tv;
        tv.tv_sec = ms / 1000;
        tv.tv_usec = (ms % 1000) * 1000;
        setsockopt(fd, SOL_SOCKET, SO_RCVTIMEO, &tv, sizeof(tv));
        setsockopt(fd, SOL_SOCKET, SO_SNDTIMEO, &tv, sizeof(tv));
      }

      bool send_message(int fd, const std::string &m, std::string *error) {
        unsigned char hdr[4];
        encode_length(static_cast<uint32_t>(m.size()), hdr);
        if (!write_all(fd, reinterpret_cast<const char *>(hdr), 4) ||
            !write_all(fd, m.data(), m.size())) {
          if (error) *error = errno_text();
          return false;
        }
        return true;
      }

      bool receive_message(int fd, std::string &m, std::string *error) {
        unsigned char hdr[4];
        if (!read_all(fd, reinterpret_cast<char *>(hdr), 4)) {
          if (error) *error = errno_text();
          return false;
        }
        const uint32_t n = decode_length(hdr);
        if (n > kMaxMessage) {
          if (error) *error = "message too large";
          return false;
        }
        m.assign(n, '\0');
        if (n && !read_all(fd, &m[0], n)) {
          if (error) *error = errno_text();
          return false;
        }
        return true;
      }

      bool fill_address(const std::string &path, sockaddr_un &addr, std::string *error) {
        std::memset(&addr, 0, sizeof(addr));
        addr.sun_family = AF_UNIX;
        if (path.size() >= sizeof(addr.sun_path)) {
          if (error) *error = "socket path too long: " + path;
          return false;
        }
        std::memcpy(addr.sun_path, path.c_str(), path.size());
        return true;
      }

      int connect_to(const std::string &path) {
        sockaddr_un addr;
        if (!fill_address(path, addr, nullptr)) return -1;
        const int fd = socket(AF_UNIX, SOCK_STREAM, 0);
        if (fd < 0) return -1;
        if (::connect(fd, reinterpret_cast<sockaddr *>(&addr), sizeof(addr)) == 0) return fd;
        ::close(fd);
        return -1;
      }
    } // namespace

    std::string default_endpoint() {
      const char *run = std::getenv("XDG_RUNTIME_DIR");
      if (run && *run) return std::string(run) + "/faster-astap-index.sock";
      const char *tmp = std::getenv("TMPDIR");
      const std::string dir = tmp && *tmp ? tmp : "/tmp";
      return dir + "/faster-astap-index-" + std::to_string(getuid()) + ".sock";
    }

    bool request(const std::string &endpoint, const std::string &message, std::string &reply,
                 int connect_timeout_ms, int reply_timeout_ms, std::string *error) {
      const auto deadline =
          std::chrono::steady_clock::now() + std::chrono::milliseconds(connect_timeout_ms);
      int fd = -1;
      for (;;) {
        fd = connect_to(endpoint);
        if (fd >= 0) break;
        if (std::chrono::steady_clock::now() >= deadline) {
          if (error) *error = "no server listening on " + endpoint;
          return false;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
      }
      set_timeout(fd, reply_timeout_ms);
      const bool ok =
          send_message(fd, message, error) && receive_message(fd, reply, error);
      ::close(fd);
      return ok;
    }

    bool endpoint_alive(const std::string &endpoint) {
      const int fd = connect_to(endpoint);
      if (fd < 0) return false;
      ::close(fd);
      return true;
    }

    struct Server::Impl {
      int listen_fd = -1;
      std::vector<std::thread> workers;
      Handler handler;
      std::atomic<bool> stopping{false};
      std::mutex done_mutex;
      std::condition_variable done;
    };

    Server::Server() : impl_(new Impl) {}

    Server::~Server() { stop(); }

    bool Server::start(const std::string &endpoint, int instances, Handler handler,
                       std::string *error) {
      endpoint_ = endpoint;
      impl_->handler = std::move(handler);
      if (instances < 1) instances = 1;

      // A socket file outlives the process that made it, so one left behind by a
      // crash would block every future server. Connecting to it is the only way
      // to tell a live server from a dead file.
      if (endpoint_alive(endpoint)) {
        if (error) *error = "another server already holds " + endpoint;
        return false;
      }
      ::unlink(endpoint.c_str());

      sockaddr_un addr;
      if (!fill_address(endpoint, addr, error)) return false;
      impl_->listen_fd = socket(AF_UNIX, SOCK_STREAM, 0);
      if (impl_->listen_fd < 0) {
        if (error) *error = "socket: " + errno_text();
        return false;
      }
      if (::bind(impl_->listen_fd, reinterpret_cast<sockaddr *>(&addr), sizeof(addr)) != 0 ||
          ::listen(impl_->listen_fd, 16) != 0) {
        if (error) *error = "bind/listen: " + errno_text();
        ::close(impl_->listen_fd);
        impl_->listen_fd = -1;
        return false;
      }
      // The solver answers whoever asks, so the socket is the access control:
      // owner only.
      ::chmod(endpoint.c_str(), S_IRUSR | S_IWUSR);
      (void) instances;
      return true;
    }

    void Server::run() {
      const int workers = 4;
      for (int i = 0; i < workers; i++) {
        impl_->workers.emplace_back([this] {
          while (!impl_->stopping) {
            const int fd = ::accept(impl_->listen_fd, nullptr, nullptr);
            if (fd < 0) continue;
            set_timeout(fd, 300000);
            std::string req;
            if (receive_message(fd, req, nullptr)) {
              const std::string rep = impl_->handler(req);
              send_message(fd, rep, nullptr);
            }
            ::close(fd);
          }
        });
      }
      std::unique_lock<std::mutex> lk(impl_->done_mutex);
      impl_->done.wait(lk, [this] { return impl_->stopping.load(); });
      lk.unlock();

      // Closing the listening socket is what makes the parked accept() calls
      // return, which is how the workers notice they are finished.
      if (impl_->listen_fd >= 0) {
        ::shutdown(impl_->listen_fd, SHUT_RDWR);
        ::close(impl_->listen_fd);
        impl_->listen_fd = -1;
      }
      for (std::thread &t : impl_->workers)
        if (t.joinable()) t.join();
      impl_->workers.clear();
      ::unlink(endpoint_.c_str());
    }

    void Server::stop() {
      if (impl_->stopping.exchange(true)) return;
      std::lock_guard<std::mutex> lk(impl_->done_mutex);
      impl_->done.notify_all();
    }
#endif
  } // namespace ipc
} // namespace astap
