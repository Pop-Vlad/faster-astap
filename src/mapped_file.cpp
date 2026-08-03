#include "astap/mapped_file.h"

#include <filesystem>

#ifdef _WIN32
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace astap {
  namespace {
    bool fail(std::string *error, const std::string &why) {
      if (error) *error = why;
      return false;
    }
  } // namespace

  MappedFile &MappedFile::operator=(MappedFile &&o) noexcept {
    if (this != &o) {
      close();
      data_ = o.data_;
      size_ = o.size_;
      file_ = o.file_;
      mapping_ = o.mapping_;
      fd_ = o.fd_;
      o.data_ = nullptr;
      o.size_ = 0;
      o.file_ = nullptr;
      o.mapping_ = nullptr;
      o.fd_ = -1;
    }
    return *this;
  }

#ifdef _WIN32

  bool MappedFile::open(const std::string &path, std::string *error) {
    close();
    // Through std::filesystem::path so a path outside the system codepage still
    // resolves; the rest of the port opens files by narrow name, which does not.
    const std::wstring wide = std::filesystem::path(path).wstring();

    // FILE_SHARE_READ only: another process must not be able to rewrite the
    // file while these pages are live, because the mapping would follow the
    // change under the reader's feet. A writer is refused instead, which the
    // cache writer already reports as a warning it can carry on from.
    HANDLE file = CreateFileW(wide.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr, OPEN_EXISTING,
                              FILE_ATTRIBUTE_NORMAL, nullptr);
    if (file == INVALID_HANDLE_VALUE) return fail(error, "cannot open " + path);

    LARGE_INTEGER size;
    if (!GetFileSizeEx(file, &size) || size.QuadPart <= 0) {
      CloseHandle(file);
      return fail(error, path + " is empty or its size cannot be read");
    }
    // A 32 bit build cannot map a file larger than its address space. Saying so
    // beats a mapping failure with no explanation.
    if (static_cast<unsigned long long>(size.QuadPart) > SIZE_MAX) {
      CloseHandle(file);
      return fail(error, path + " is too large to map in a 32 bit process");
    }

    HANDLE mapping = CreateFileMappingW(file, nullptr, PAGE_READONLY, 0, 0, nullptr);
    if (mapping == nullptr) {
      CloseHandle(file);
      return fail(error, "cannot map " + path);
    }
    void *view = MapViewOfFile(mapping, FILE_MAP_READ, 0, 0, 0);
    if (view == nullptr) {
      CloseHandle(mapping);
      CloseHandle(file);
      return fail(error, "cannot map " + path);
    }

    // Both handles are dropped here: the section keeps the file open for as long
    // as a view of it exists, so the view outlives them. Holding no descriptor
    // per mapping is what lets a caller keep every area of a 1476 file database
    // mapped at once without running into a descriptor limit.
    CloseHandle(mapping);
    CloseHandle(file);

    data_ = static_cast<const uint8_t *>(view);
    size_ = static_cast<size_t>(size.QuadPart);
    return true;
  }

  void MappedFile::close() {
    if (data_) UnmapViewOfFile(static_cast<LPCVOID>(data_));
    data_ = nullptr;
    size_ = 0;
    file_ = nullptr;
    mapping_ = nullptr;
  }

  // Windows offers PrefetchVirtualMemory, which is the opposite request, and no
  // way to say "do not read ahead". Nothing to do.
  void MappedFile::advise_random() const {}

#else

  bool MappedFile::open(const std::string &path, std::string *error) {
    close();
    const int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) return fail(error, "cannot open " + path);

    struct stat st;
    if (::fstat(fd, &st) != 0 || st.st_size <= 0) {
      ::close(fd);
      return fail(error, path + " is empty or its size cannot be read");
    }
    if (static_cast<unsigned long long>(st.st_size) > SIZE_MAX) {
      ::close(fd);
      return fail(error, path + " is too large to map in a 32 bit process");
    }

    void *view = ::mmap(nullptr, static_cast<size_t>(st.st_size), PROT_READ, MAP_SHARED, fd, 0);
    if (view == MAP_FAILED) {
      ::close(fd);
      return fail(error, "cannot map " + path);
    }

    // The descriptor is dropped here: a mapping keeps its own reference to the
    // file, so the pages outlive it. Holding no descriptor per mapping is what
    // lets a caller keep every area of a 1476 file database mapped at once
    // without running into RLIMIT_NOFILE, which defaults to 1024.
    ::close(fd);

    data_ = static_cast<const uint8_t *>(view);
    size_ = static_cast<size_t>(st.st_size);
    return true;
  }

  void MappedFile::close() {
    if (data_) ::munmap(const_cast<uint8_t *>(data_), size_);
    data_ = nullptr;
    size_ = 0;
    fd_ = -1;
  }

  void MappedFile::advise_random() const {
#ifdef MADV_RANDOM
    if (data_) ::madvise(const_cast<uint8_t *>(data_), size_, MADV_RANDOM);
#endif
  }

#endif
} // namespace astap
