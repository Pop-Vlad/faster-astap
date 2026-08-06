// A read-only memory mapping of a whole file.
//
// Both large files this solver reads — the index cache and the star database —
// are touched in scattered pieces rather than end to end, which is what a
// mapping is for. Against reading into a buffer it saves three things: the copy,
// the resident memory that copy occupies, and, because the page cache is what
// backs the mapping, a second copy when another process opens the same file.
// The pages a solve never touches are never moved at all.
//
// Read only, and the whole file at once: neither caller writes through a mapping
// or needs a window into one, and leaving both out keeps each platform half to a
// dozen lines.

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>

namespace astap {
  class MappedFile {
  public:
    MappedFile() = default;
    ~MappedFile() { close(); }

    // Movable, not copyable: a mapping is a handle, and two objects unmapping
    // the same address is a double free. Share one through a shared_ptr.
    MappedFile(const MappedFile &) = delete;
    MappedFile &operator=(const MappedFile &) = delete;
    MappedFile(MappedFile &&o) noexcept { *this = std::move(o); }
    MappedFile &operator=(MappedFile &&o) noexcept;

    // Maps `path` in its entirety. An empty file is reported as a failure:
    // there would be nothing to point at, and every caller wants a pointer.
    bool open(const std::string &path, std::string *error = nullptr);
    void close();

    bool is_open() const { return data_ != nullptr; }
    const uint8_t *data() const { return data_; }
    size_t size() const { return size_; }

  private:
    const uint8_t *data_ = nullptr;
    size_t size_ = 0;
    // HANDLEs on Windows, held as void * so this header does not drag windows.h
    // into everything that includes it. The POSIX half needs only the descriptor.
    void *file_ = nullptr;
    void *mapping_ = nullptr;
    int fd_ = -1;
  };
} // namespace astap
