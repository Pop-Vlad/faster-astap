#include "astap/parallel.h"

#include <condition_variable>
#include <mutex>

#ifdef __linux__
#include <sched.h>
#endif

namespace astap {
  namespace {
    unsigned g_threads = 0;

    // Number of hardware threads the process is actually allowed to run on.
    //
    // This counts logical processors, SMT siblings included: on a blind solve the
    // spiral search keeps every sibling busy with independent work, and using all
    // 20 threads of a 14 core machine measured about 14% faster than using 14.
    //
    // hardware_concurrency() reports the size of the machine rather than what this
    // process may use, so the affinity mask is preferred where available. That is
    // what makes the count correct under taskset, inside a container, or on a
    // batch scheduler that pins jobs to a subset of the CPUs.
    unsigned available_threads() {
#ifdef __linux__
      cpu_set_t mask;
      CPU_ZERO(&mask);
      if (sched_getaffinity(0, sizeof(mask), &mask) == 0) {
        const int n = CPU_COUNT(&mask);
        if (n > 0) return static_cast<unsigned>(n);
      }
#endif
      const unsigned hw = std::thread::hardware_concurrency();
      return hw == 0 ? 1 : hw;
    }

    // A persistent pool. The spiral search runs thousands of small parallel
    // batches, and creating a thread per batch costs more than the work itself.
    class Pool {
    public:
      static Pool &instance() {
        static Pool p;
        return p;
      }

      unsigned size() const { return workers_ + 1; } // workers plus the calling thread

      // Runs body(chunk_lo, chunk_hi, slot) for `n` contiguous chunks. Blocks until
      // all chunks are done. Not reentrant: nested calls run sequentially.
      void run(const std::function<void(unsigned)> &body, unsigned n) {
        if (n <= 1 || busy_) {
          for (unsigned i = 0; i < n; i++) body(i);
          return;
        }
        ensure_started();

        {
          std::lock_guard<std::mutex> lk(m_);
          body_ = &body;
          chunks_ = n;
          next_ = 1; // chunk 0 runs on the calling thread
          done_ = 0;
          busy_ = true;
            }
        cv_work_.notify_all();

        body(0); // the caller takes the first chunk

        // Help with whatever is left, then wait for the stragglers.
        for (;;) {
          unsigned mine;
          {
            std::lock_guard<std::mutex> lk(m_);
            if (next_ >= chunks_) break;
            mine = next_++;
          }
          body(mine);
          std::lock_guard<std::mutex> lk(m_);
          done_++;
        }

        std::unique_lock<std::mutex> lk(m_);
        cv_done_.wait(lk, [this] { return done_ >= chunks_ - 1; });
        busy_ = false;
        body_ = nullptr;
      }

      ~Pool() {
        {
          std::lock_guard<std::mutex> lk(m_);
          stop_ = true;
            }
        cv_work_.notify_all();
        for (std::thread &t: pool_) t.join();
      }

    private:
      void ensure_started() {
        const unsigned want = thread_count() > 0 ? thread_count() - 1 : 0;
        if (workers_ == want) return;
        // Workers are only ever added; a lower thread count simply leaves some idle.
        for (unsigned i = workers_; i < want; i++)
          pool_.emplace_back([this] { worker_loop(); });
        workers_ = want;
      }

      void worker_loop() {
      for (;;) {
        unsigned mine;
        const std::function<void(unsigned)> *body;
        {
          std::unique_lock<std::mutex> lk(m_);
          cv_work_.wait(lk, [this] { return stop_ || (busy_ && next_ < chunks_); });
          if (stop_) return;
          // The task pointer is read under the same lock that hands out the
          // chunk. Caching it across iterations would be a use after return:
          // it points at the caller's stack, and the caller may already have
          // returned and started the next batch.
          mine = next_++;
          body = body_;
        }
        (*body)(mine);
        {
          std::lock_guard<std::mutex> lk(m_);
          done_++;
        }
        cv_done_.notify_one();
      }
    }

    std::vector<std::thread> pool_;
      unsigned workers_ = 0;
      std::mutex m_;
      std::condition_variable cv_work_, cv_done_;
      const std::function<void(unsigned)> *body_ = nullptr;
      unsigned chunks_ = 0, next_ = 0, done_ = 0;
        bool busy_ = false, stop_ = false;
    };
  } // namespace

  void set_thread_count(unsigned n) { g_threads = n; }

  unsigned thread_count() {
    if (g_threads != 0) return g_threads;
    static const unsigned n = available_threads();
    return n == 0 ? 1 : n;
  }

  unsigned range_chunks(size_t count) {
    if (count == 0) return 0;
    return static_cast<unsigned>(std::min<size_t>(thread_count(), count));
  }

  void parallel_ranges(size_t begin, size_t end,
                       const std::function<void(size_t, size_t, unsigned)> &body) {
    if (end <= begin) return;
    const size_t count = end - begin;
    const unsigned n = range_chunks(count);
    if (n <= 1) {
      body(begin, end, 0);
      return;
    }

    const size_t chunk = count / n;
    const size_t rest = count % n;

    // The first `rest` chunks get one extra item, keeping the chunks contiguous
    // and in increasing order.
    auto run_chunk = [&](unsigned t) {
      const size_t lo = begin + t * chunk + std::min<size_t>(t, rest);
      const size_t hi = lo + chunk + (t < rest ? 1 : 0);
      body(lo, hi, t);
    };

    Pool::instance().run(run_chunk, n);
  }
} // namespace astap
