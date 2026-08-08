// A local request/response channel between the resident solver and the small
// client that an imaging application launches.
//
// Deliberately the least machinery that does the job: one connection carries one
// request and one response, each a length-prefixed UTF-8 blob, and the transport
// is whatever the platform already offers for talking to another process on the
// same machine — a named pipe on Windows, a Unix domain socket elsewhere. No
// sockets on a TCP port, so nothing is reachable from off the machine and no
// firewall prompt appears the first time an observatory PC runs it.
//
// The payload format is not the transport's business; see index_server.cpp for
// the key=value messages that actually travel over it.

#pragma once

#include <atomic>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace astap {
  namespace ipc {
    // Where a server listens and a client looks, when neither was told
    // otherwise. `\\.\pipe\faster-astap-index` on Windows, a socket under the
    // user's runtime directory elsewhere.
    std::string default_endpoint();

    // Sends one request and waits for the reply. `connect_timeout_ms` covers a
    // server that is still starting up, `reply_timeout_ms` a solve that is
    // taking its time; both are bounded so that a wedged server costs a client
    // an error rather than a night.
    //
    // Returns false and fills `error` when the server could not be reached,
    // which is the case the caller has to tell apart from a refused request: the
    // first is worth starting a server for, the second is not.
    bool request(const std::string &endpoint, const std::string &message, std::string &reply,
                 int connect_timeout_ms, int reply_timeout_ms, std::string *error);

    // True when something is listening on `endpoint` right now. Cheaper than a
    // request and used to decide whether a server has to be started.
    bool endpoint_alive(const std::string &endpoint);

    // Serves requests until `stop` is called.
    //
    // `handler` runs on one of `instances` worker threads, so a slow request
    // does not block the next client from connecting. It says nothing about
    // whether the work itself may overlap — the solver cannot, and serialises
    // itself with a mutex.
    class Server {
    public:
      using Handler = std::function<std::string(const std::string & request)>;

      Server();

      ~Server();

      Server(const Server &) = delete;

      Server &operator=(const Server &) = delete;

      // Takes the endpoint. Fails when another server already holds it, which is
      // what keeps two servers from half-answering each other's clients.
      bool start(const std::string &endpoint, int instances, Handler handler, std::string *error);

      // Blocks until stop() is called from another thread or from a handler.
      void run();

      void stop();

      const std::string &endpoint() const { return endpoint_; }

    private:
      struct Impl;
      std::unique_ptr<Impl> impl_;
      std::string endpoint_;
    };
  } // namespace ipc
} // namespace astap
