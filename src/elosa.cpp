// elosa: Exhaustive Lesion Overlap Space Analysis.
//
// Enumerates every "group": a set of patients G such that the region
// R(G) = {voxels lesioned in every patient of G} is non-empty and no other
// patient's lesion covers all of R(G). Each group corresponds to exactly one
// distinct overlap region, so the number of groups is the number of distinct
// regions a lesion dataset can distinguish.
//
// This is closed-itemset mining: voxels with the same pattern of lesioned
// patients are merged into weighted "atoms" (transactions), patients are
// items, and groups are the closed itemsets. The search is LCM-style
// (Uno et al. 2004): prefix-preserving closure extension, so each group is
// generated exactly once with no duplicate checks; occurrence lists built by
// counting sort; per-node database reduction (items that can no longer be
// added are dropped and identical atoms merged); bit-parallel closure tests
// with early exit. Work is shared between threads by handing sub-trees to
// idle threads.
//
// Build: see README.md (g++ -O3 -march=native -std=c++17 -pthread).
#define _CRT_SECURE_NO_WARNINGS   // MSVC: allow fopen
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>

#if defined(_MSC_VER)
#include <intrin.h>
static inline int popcount64(uint64_t x) { return (int)__popcnt64(x); }
static inline int ctz64(uint64_t x) { unsigned long i; _BitScanForward64(&i, x); return (int)i; }
static inline int clz64(uint64_t x) { unsigned long i; _BitScanReverse64(&i, x); return 63 - (int)i; }
#else
static inline int popcount64(uint64_t x) { return __builtin_popcountll(x); }
static inline int ctz64(uint64_t x) { return __builtin_ctzll(x); }
static inline int clz64(uint64_t x) { return __builtin_clzll(x); }
#endif

// ------------------------------------------------------------------ threads
// std::thread/std::mutex are missing from MinGW builds that use the "win32"
// thread model (e.g. some MATLAB MinGW-w64 installs), so on those use a small
// Win32 implementation with the same interface. Define ELOSA_WIN32_THREADS to
// force it, or ELOSA_STD_THREADS to force the standard library.
#if !defined(ELOSA_STD_THREADS) && !defined(ELOSA_WIN32_THREADS) && defined(_WIN32) && \
    defined(__GLIBCXX__) && !defined(_GLIBCXX_HAS_GTHREADS)
#define ELOSA_WIN32_THREADS
#endif

#ifdef ELOSA_WIN32_THREADS
#include <functional>
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
namespace esync {
class mutex {
public:
    mutex() { InitializeSRWLock(&l_); }
    mutex(const mutex&) = delete;
    mutex& operator=(const mutex&) = delete;
    void lock() { AcquireSRWLockExclusive(&l_); }
    void unlock() { ReleaseSRWLockExclusive(&l_); }
    SRWLOCK* native() { return &l_; }
private:
    SRWLOCK l_;
};
template <class M> class lock_guard {
public:
    explicit lock_guard(M& m) : m_(m) { m_.lock(); }
    ~lock_guard() { m_.unlock(); }
    lock_guard(const lock_guard&) = delete;
    lock_guard& operator=(const lock_guard&) = delete;
private:
    M& m_;
};
template <class M> class unique_lock {
public:
    explicit unique_lock(M& m) : m_(&m) { m_->lock(); }
    ~unique_lock() { m_->unlock(); }
    unique_lock(const unique_lock&) = delete;
    unique_lock& operator=(const unique_lock&) = delete;
    M* mutex() { return m_; }
private:
    M* m_;
};
class condition_variable {
public:
    condition_variable() { InitializeConditionVariable(&c_); }
    condition_variable(const condition_variable&) = delete;
    condition_variable& operator=(const condition_variable&) = delete;
    void notify_one() { WakeConditionVariable(&c_); }
    void notify_all() { WakeAllConditionVariable(&c_); }
    template <class Pred> void wait(unique_lock<esync::mutex>& lk, Pred pred) {
        while (!pred()) SleepConditionVariableSRW(&c_, lk.mutex()->native(), INFINITE, 0);
    }
    // returns pred() (true = condition met, false = timed out)
    template <class Rep, class Period, class Pred>
    bool wait_for(unique_lock<esync::mutex>& lk, const std::chrono::duration<Rep, Period>& d, Pred pred) {
        const auto end = std::chrono::steady_clock::now() + d;
        while (!pred()) {
            const auto now = std::chrono::steady_clock::now();
            if (now >= end) return pred();
            const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - now).count();
            SleepConditionVariableSRW(&c_, lk.mutex()->native(), DWORD(ms > 0 ? ms : 1), 0);
        }
        return true;
    }
private:
    CONDITION_VARIABLE c_;
};
class thread {
public:
    template <class F> explicit thread(F f) {
        auto* fn = new std::function<void()>(std::move(f));
        h_ = CreateThread(nullptr, 0, &thread::entry, fn, 0, nullptr);
        if (!h_) { delete fn; throw std::runtime_error("cannot create thread"); }
    }
    thread(thread&& o) noexcept : h_(o.h_) { o.h_ = nullptr; }
    thread& operator=(thread&& o) noexcept { std::swap(h_, o.h_); return *this; }
    thread(const thread&) = delete;
    ~thread() { if (h_) CloseHandle(h_); }
    void join() { WaitForSingleObject(h_, INFINITE); CloseHandle(h_); h_ = nullptr; }
    static unsigned hardware_concurrency() {
        SYSTEM_INFO si;
        GetSystemInfo(&si);
        return unsigned(si.dwNumberOfProcessors);
    }
private:
    static DWORD WINAPI entry(LPVOID p) {
        std::unique_ptr<std::function<void()>> fn(static_cast<std::function<void()>*>(p));
        (*fn)();
        return 0;
    }
    HANDLE h_ = nullptr;
};
}  // namespace esync
#else
#include <condition_variable>
#include <mutex>
#include <thread>
namespace esync {
using std::condition_variable;
using std::lock_guard;
using std::mutex;
using std::thread;
using std::unique_lock;
}  // namespace esync
#endif

// Is the process with this id still running? Used so that elosa stops if
// the MATLAB session that started it is killed or crashes.
#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
static bool process_alive(unsigned long pid) {
    HANDLE h = OpenProcess(SYNCHRONIZE, FALSE, DWORD(pid));
    if (!h) return GetLastError() == ERROR_ACCESS_DENIED;   // exists but not ours to open
    const DWORD r = WaitForSingleObject(h, 0);
    CloseHandle(h);
    return r == WAIT_TIMEOUT;
}
#else
#include <cerrno>
#include <signal.h>
static bool process_alive(unsigned long pid) {
    return kill(pid_t(pid), 0) == 0 || errno == EPERM;
}
#endif

namespace {

// ------------------------------------------------------------------ options
struct Options {
    std::string input;
    std::string out;          // binary file of groups ("" = none)
    std::string summary;      // summary text file ("" = stdout only)
    int threads = 0;          // 0 = hardware concurrency
    uint64_t min_voxels = 1;  // minimum region size
    int max_patients = 0;     // 0 = no limit
    bool quiet = false;
    bool covering = false;
    std::string progress;     // file to rewrite with progress every second ("" = none)
    unsigned long parent_pid = 0;   // stop if this process ends (0 = don't watch)    // also list non-grouped patients covering each region
    std::string order = "ascending";  // item order: ascending / descending / input
    long estimate = 0;        // > 0: estimate the number of groups from this many random walks
    uint64_t seed = 1;
};

[[noreturn]] void usage(const char* msg = nullptr) {
    if (msg) std::cerr << "elosa: " << msg << "\n\n";
    std::cerr <<
        "usage: elosa INPUT [options]\n"
        "  INPUT                 binary lesion file written by elosa_write_input.m\n"
        "  --out PREFIX          write every group to PREFIX.groups and PREFIX.members\n"
        "                        (read with elosa_read_groups.m)\n"
        "  --covering            with --out, also list the non-grouped patients whose\n"
        "                        lesions cover each region (for PPV on unimpaired/test patients)\n"
        "  --summary FILE        write counts and histograms to FILE\n"
        "  --threads N           worker threads (default: all cores)\n"
        "  --min-voxels N        only regions of at least N voxels (default 1 = exhaustive)\n"
        "  --max-patients N      only groups of at most N patients (default: no limit)\n"
        "  --order ascending|descending|input   patient order used by the search\n"
        "                        (affects speed only, never the result; default ascending)\n"
        "  --estimate N          don't enumerate: estimate the number of groups from N random\n"
        "                        walks down the search tree (Knuth's estimator)\n"
        "  --seed S              random seed for --estimate (default 1)\n"
        "  --progress FILE       rewrite FILE every second with: phase, fraction done,\n"
        "                        groups so far, seconds elapsed, rough seconds left\n"
        "  --parent-pid PID      stop as soon as process PID ends (elosa_run passes MATLAB's)\n"
        "  --quiet               no progress messages\n";
    std::exit(2);
}

Options parse_args(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) usage(("missing value for " + a).c_str());
            return argv[++i];
        };
        if (a == "--out") o.out = next();
        else if (a == "--summary") o.summary = next();
        else if (a == "--threads") o.threads = std::stoi(next());
        else if (a == "--min-voxels") o.min_voxels = std::stoull(next());
        else if (a == "--max-patients") o.max_patients = std::stoi(next());
        else if (a == "--order") o.order = next();
        else if (a == "--quiet") o.quiet = true;
        else if (a == "--covering") o.covering = true;
        else if (a == "--progress") o.progress = next();
        else if (a == "--parent-pid") o.parent_pid = std::stoul(next());
        else if (a == "--estimate") o.estimate = std::stol(next());
        else if (a == "--seed") o.seed = std::stoull(next());
        else if (a == "-h" || a == "--help") usage();
        else if (!a.empty() && a[0] == '-') usage(("unknown option " + a).c_str());
        else if (o.input.empty()) o.input = a;
        else usage("more than one input file");
    }
    if (o.input.empty()) usage("no input file");
    if (o.order != "ascending" && o.order != "descending" && o.order != "input")
        usage("--order must be ascending, descending or input");
    if (o.min_voxels < 1) o.min_voxels = 1;
    return o;
}

// ------------------------------------------------------------------ input
// File layout (little-endian):
//   char[8]  "ELOSAIN1"
//   uint32   n (patients), uint32 V (voxels)
//   uint8    use[n]    1 = patient may be in groups (e.g. impaired), 0 = only
//                      counted when its lesion covers a region
//   uint8    lesion[V][n]   voxel-major (MATLAB column-major n x V), 0/1
struct Input {
    uint32_t n = 0, V = 0;
    std::vector<uint8_t> use;
    std::vector<uint8_t> lesion;  // V*n
};

Input read_input(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path);
    char magic[8];
    f.read(magic, 8);
    if (!f || std::memcmp(magic, "ELOSAIN1", 8) != 0) throw std::runtime_error("not an ELOSA input file: " + path);
    Input in;
    f.read(reinterpret_cast<char*>(&in.n), 4);
    f.read(reinterpret_cast<char*>(&in.V), 4);
    if (!f || in.n == 0) throw std::runtime_error("bad header in " + path);
    in.use.resize(in.n);
    f.read(reinterpret_cast<char*>(in.use.data()), in.n);
    in.lesion.resize(size_t(in.n) * in.V);
    f.read(reinterpret_cast<char*>(in.lesion.data()), std::streamsize(in.lesion.size()));
    if (!f) throw std::runtime_error("file too short: " + path);
    return in;
}

// ------------------------------------------------------------------ database
// A (conditional) database: T weighted atoms, each with
//   bits  - the atom's full patient bitset (AND of merged atoms), all patients
//   items - the patients that may still be added (sorted, > current core)
// Buffers only ever grow (T and NI are the sizes in use), so the search does
// not allocate or zero memory once it has warmed up.
struct DB {
    int W = 0;
    size_t T = 0, NI = 0;
    std::vector<uint64_t> bits;
    std::vector<uint64_t> weight;
    std::vector<uint32_t> off{0};
    std::vector<uint32_t> items;
    std::vector<uint64_t> hash;     // per atom, hash of its item list (filled by the builder)
    size_t size() const { return T; }
    const uint64_t* row(size_t t) const { return bits.data() + t * W; }
    uint64_t* row(size_t t) { return bits.data() + t * W; }
    void reserve(size_t t, size_t ni) {
        if (bits.size() < t * W) bits.resize(t * W);
        if (weight.size() < t) weight.resize(t);
        if (off.size() < t + 1) off.resize(t + 1);
        if (items.size() < ni) items.resize(ni);
        if (hash.size() < t) hash.resize(t);
        off[0] = 0;
    }
    // copy of just the part in use (for handing work to another thread)
    DB compact() const {
        DB d;
        d.W = W; d.T = T; d.NI = NI;
        d.bits.assign(bits.begin(), bits.begin() + T * W);
        d.weight.assign(weight.begin(), weight.begin() + T);
        d.off.assign(off.begin(), off.begin() + T + 1);
        d.items.assign(items.begin(), items.begin() + NI);
        return d;
    }
};

static inline uint64_t hash_step(uint64_t h, uint32_t x) { return (h ^ x) * 0xff51afd7ed558ccdULL; }
static inline uint64_t hash_final(uint64_t h) { return h ^ (h >> 29) ^ (h >> 47); }

struct MergeTable {
    std::vector<uint32_t> slot, stamp;
    uint32_t gen = 0;
};

// Merge atoms whose item lists are identical (weights add, bitsets AND), in
// place and in one pass with an open-addressing hash table.
void merge_duplicates(DB& db, MergeTable& h) {
    const size_t T = db.T;
    if (T < 2) return;
    size_t cap = 16;
    while (cap < 2 * T) cap <<= 1;
    if (h.slot.size() < cap) { h.slot.assign(cap, 0); h.stamp.assign(cap, 0); h.gen = 0; }
    if (++h.gen == 0) { std::fill(h.stamp.begin(), h.stamp.end(), 0); h.gen = 1; }
    const size_t mask = cap - 1;
    const int W = db.W;
    uint32_t* items = db.items.data();
    uint32_t* off = db.off.data();
    size_t nt = 0;
    uint32_t start = off[0];
    for (size_t t = 0; t < T; ++t) {
        const uint32_t end = off[t + 1];          // read before off[nt + 1] may overwrite it
        const uint32_t len = end - start;
        size_t pos = hash_final(db.hash[t]) & mask;
        for (;;) {
            if (h.stamp[pos] != h.gen) {          // new list: move it down to slot nt
                h.stamp[pos] = h.gen;
                h.slot[pos] = uint32_t(nt);
                const uint32_t dst = off[nt];
                if (nt != t) {
                    std::memmove(items + dst, items + start, len * sizeof(uint32_t));
                    std::memcpy(db.row(nt), db.row(t), W * sizeof(uint64_t));
                    db.weight[nt] = db.weight[t];
                }
                off[nt + 1] = dst + len;
                ++nt;
                break;
            }
            const uint32_t u = h.slot[pos];
            if (off[u + 1] - off[u] == len && std::memcmp(items + off[u], items + start, len * sizeof(uint32_t)) == 0) {
                db.weight[u] += db.weight[t];
                uint64_t* r = db.row(u);
                const uint64_t* q = db.row(t);
                for (int w = 0; w < W; ++w) r[w] &= q[w];
                break;
            }
            pos = (pos + 1) & mask;
        }
        start = end;
    }
    db.T = nt;
    db.NI = off[nt];
}

// ------------------------------------------------------------------ results
struct Stats {
    uint64_t groups = 0;
    std::vector<uint64_t> by_size;        // [k] groups with k patients
    std::vector<uint64_t> by_region_log2; // [b] groups with region size in [2^b, 2^(b+1))
    void init(int m) { by_size.assign(m + 1, 0); by_region_log2.assign(64, 0); }
    void add(const Stats& o) {
        groups += o.groups;
        for (size_t i = 0; i < by_size.size(); ++i) by_size[i] += o.by_size[i];
        for (size_t i = 0; i < 64; ++i) by_region_log2[i] += o.by_region_log2[i];
    }
};

// Output: PREFIX.groups holds one fixed-size record per group,
//   uint32 [patients in group, other covering patients, region voxels lo, hi]
// after a header (char[8] "ELOSAGR2", uint32 n, uint32 flags); PREFIX.members
// holds the 1-based patient indices of each group in the same order, followed
// (if flags & 1, --covering) by the other patients whose lesions cover the region.
// Indices are uint16 if flags & 2 (fewer than 65536 patients), else uint32.
struct OutBuf {
    std::vector<uint32_t> head, mem;
};

class GroupWriter {
public:
    GroupWriter(const std::string& prefix, uint32_t n, uint32_t flags) : narrow_((flags & 2) != 0) {
        g_ = std::fopen((prefix + ".groups").c_str(), "wb");
        m_ = std::fopen((prefix + ".members").c_str(), "wb");
        if (!g_ || !m_) throw std::runtime_error("cannot write " + prefix + ".groups/.members");
        std::fwrite("ELOSAGR2", 1, 8, g_);
        std::fwrite(&n, 4, 1, g_);
        std::fwrite(&flags, 4, 1, g_);
    }
    ~GroupWriter() {
        if (g_) std::fclose(g_);
        if (m_) std::fclose(m_);
    }
    void write(OutBuf& b) {
        esync::lock_guard<esync::mutex> lock(mu_);
        if (std::fwrite(b.head.data(), 4, b.head.size(), g_) != b.head.size() ||
            !write_members(b.mem))
            failed_ = true;
        b.head.clear();
        b.mem.clear();
    }
    bool failed() const { return failed_; }
private:
    bool write_members(const std::vector<uint32_t>& mem) {
        if (!narrow_) return std::fwrite(mem.data(), 4, mem.size(), m_) == mem.size();
        narrow_buf_.assign(mem.begin(), mem.end());
        return std::fwrite(narrow_buf_.data(), 2, narrow_buf_.size(), m_) == narrow_buf_.size();
    }
    bool narrow_;
    std::vector<uint16_t> narrow_buf_;
    std::FILE* g_ = nullptr;
    std::FILE* m_ = nullptr;
    esync::mutex mu_;
    bool failed_ = false;
};

// ------------------------------------------------------------------ search
struct Task {
    std::vector<uint64_t> G;  // closed set (bitset over all patients; only "use" bits set)
    int core;
    DB db;
    double share = 1.0;       // this sub-tree's share of the whole search (for progress)
};

class Search {
public:
    Search(const Options& o, int m, int W, std::vector<uint32_t> original_index, GroupWriter* writer, int threads)
        : opt_(o), m_(m), W_(W), orig_(std::move(original_index)), writer_(writer) {
        use_mask_.assign(W, 0);
        for (int i = 0; i < m; ++i) use_mask_[i >> 6] |= 1ULL << (i & 63);
        for (int t = 0; t < threads; ++t) done_.push_back(std::make_unique<std::atomic<double>>(0.0));
    }

    // Fraction of the search finished (0-1); see process(). A rough guide only.
    double fraction() const {
        double f = 0;
        for (const auto& d : done_) f += d->load(std::memory_order_relaxed);
        return std::min(1.0, f);
    }

    void run(Task root, int threads, Stats& total) {
        pending_ = 1;
        queue_.push_back(std::make_unique<Task>(std::move(root)));
        std::vector<esync::thread> pool;
        std::vector<Stats> stats(threads);
        for (int t = 0; t < threads; ++t) pool.emplace_back([this, t, &stats] { worker(stats[t], t); });
        for (auto& th : pool) th.join();
        total.init(m_);
        for (auto& s : stats) total.add(s);
    }

    void emit_root(const std::vector<uint64_t>& G, uint64_t sup, const uint64_t* acc, Stats& s,
                   OutBuf& buf) { emit(G.data(), sup, acc, s, buf); }

private:
    struct Cand { uint32_t item, begin, end; uint64_t sup; };
    struct Worker {
        Stats stats;
        OutBuf buf;
        std::vector<uint32_t> cnt, fill, touched;
        std::vector<uint64_t> wsum;
        // per-depth buffers; deque so growing never moves the ones in use
        struct Level {
            std::vector<uint32_t> occ, occpos;   // atom index and position of the item in its list
            std::vector<Cand> cands;
            std::vector<uint64_t> acc, G;
            DB child;
        };
        std::deque<Level> levels;
        MergeTable table;
        double done = 0;                     // share of the search finished by this thread
        std::atomic<double>* published = nullptr;
    };

    void credit(Worker& w, double share) {
        w.done += share;
        w.published->store(w.done, std::memory_order_relaxed);
    }

    void worker(Stats& out, int index) {
        Worker w;
        w.published = done_[index].get();
        w.stats.init(m_);
        w.cnt.assign(m_, 0);
        w.wsum.assign(m_, 0);
        w.fill.assign(m_, 0);
        for (;;) {
            std::unique_ptr<Task> task;
            {
                esync::unique_lock<esync::mutex> lock(qm_);
                ++idle_;
                qcv_.wait(lock, [&] { return !queue_.empty() || pending_ == 0; });
                --idle_;
                if (queue_.empty()) break;
                task = std::move(queue_.front());
                queue_.pop_front();
            }
            process(w, task->G.data(), task->core, task->db, 0, task->share);
            task.reset();
            {
                esync::lock_guard<esync::mutex> lock(qm_);
                if (--pending_ == 0) qcv_.notify_all();
            }
        }
        flush(w.buf);
        out = std::move(w.stats);
    }

    void flush(OutBuf& buf) {
        if (writer_ && !buf.head.empty()) writer_->write(buf);
    }

    // record one group: [k, region voxels (lo, hi), covering non-group patients, patient indices...]
    void emit(const uint64_t* G, uint64_t sup, const uint64_t* acc, Stats& s, OutBuf& buf) {
        int k = 0, other = 0;
        for (int w = 0; w < W_; ++w) {
            k += popcount64(G[w]);
            other += popcount64(acc[w] & ~use_mask_[w]);
        }
        if ((++s.groups & 0xFFF) == 0) progress_.fetch_add(0x1000, std::memory_order_relaxed);
        ++s.by_size[k];
        const int b = 63 - clz64(sup);   // sup >= 1
        ++s.by_region_log2[b];
        if (writer_) {
            buf.head.push_back(uint32_t(k));
            buf.head.push_back(uint32_t(other));
            buf.head.push_back(uint32_t(sup & 0xffffffffULL));
            buf.head.push_back(uint32_t(sup >> 32));
            for (int w = 0; w < W_; ++w)
                for (uint64_t x = G[w]; x; x &= x - 1)
                    buf.mem.push_back(orig_[(w << 6) + ctz64(x)]);
            if (opt_.covering)
                for (int w = 0; w < W_; ++w)
                    for (uint64_t x = acc[w] & ~use_mask_[w]; x; x &= x - 1)
                        buf.mem.push_back(orig_[(w << 6) + ctz64(x)]);
            if (buf.head.size() + buf.mem.size() > (1u << 20)) flush(buf);
        }
    }

    // Occurrence lists (counting sort) and candidate items for one database.
    void prepare(Worker& w, const DB& db, Worker::Level& L) {
        const size_t T = db.size();
        w.touched.clear();
        for (size_t t = 0; t < T; ++t) {
            for (uint32_t p = db.off[t]; p < db.off[t + 1]; ++p) {
                uint32_t i = db.items[p];
                if (w.cnt[i]++ == 0) w.touched.push_back(i);
                w.wsum[i] += db.weight[t];
            }
        }
        std::sort(w.touched.begin(), w.touched.end());
        uint32_t total = 0;
        L.cands.clear();
        for (uint32_t i : w.touched) {
            w.fill[i] = total;
            // an item is a candidate only if its region can be big enough
            if (w.wsum[i] >= opt_.min_voxels) L.cands.push_back({i, total, total + w.cnt[i], w.wsum[i]});
            total += w.cnt[i];
        }
        if (L.occ.size() < total) { L.occ.resize(total); L.occpos.resize(total); }
        for (size_t t = 0; t < T; ++t)
            for (uint32_t p = db.off[t]; p < db.off[t + 1]; ++p) {
                const uint32_t f = w.fill[db.items[p]]++;
                L.occ[f] = uint32_t(t);
                L.occpos[f] = p;
            }
        for (uint32_t i : w.touched) { w.cnt[i] = 0; w.wsum[i] = 0; }
        L.acc.resize(W_);
        L.G.resize(W_);
    }

    // Closure of G + {c.item}: the AND of every atom containing the item, into
    // L.acc and (grouped patients only) L.G. Returns the group size, or -1 if
    // the closure adds a patient j < item outside G (not a prefix-preserving
    // extension: that group is generated elsewhere) or exceeds --max-patients.
    // The words up to the item are ANDed first so rejected candidates stop early.
    int closure(const uint64_t* G, const DB& db, Worker::Level& L, const Cand& c) {
        uint64_t* acc = L.acc.data();
        const uint32_t i = c.item;
        const int wi = int(i >> 6);
        const uint64_t lowmask = (i & 63) ? ((1ULL << (i & 63)) - 1) : 0;
        const uint64_t* first = db.row(L.occ[c.begin]);
        for (int x = 0; x <= wi; ++x) acc[x] = first[x];
        for (uint32_t p = c.begin + 1; p < c.end; ++p) {
            const uint64_t* r = db.row(L.occ[p]);
            for (int x = 0; x <= wi; ++x) acc[x] &= r[x];
        }
        uint64_t bad = acc[wi] & ~G[wi] & use_mask_[wi] & lowmask;
        for (int x = 0; x < wi && !bad; ++x) bad = acc[x] & ~G[x] & use_mask_[x];
        if (bad) return -1;
        if (wi + 1 < W_) {
            for (int x = wi + 1; x < W_; ++x) acc[x] = first[x];
            for (uint32_t p = c.begin + 1; p < c.end; ++p) {
                const uint64_t* r = db.row(L.occ[p]);
                for (int x = wi + 1; x < W_; ++x) acc[x] &= r[x];
            }
        }
        int k = 0;
        for (int x = 0; x < W_; ++x) { L.G[x] = acc[x] & use_mask_[x]; k += popcount64(L.G[x]); }
        if (opt_.max_patients > 0 && k > opt_.max_patients) return -1;   // children only get bigger
        return k;
    }

    // Child database: atoms containing c.item, keeping only items > c.item that
    // are not in the new group, with identical atoms merged. Returns its size.
    size_t build_child(Worker& w, const DB& db, Worker::Level& L, const Cand& c) {
        const uint64_t* acc = L.acc.data();
        DB& child = L.child;
        child.W = W_;
        size_t cap = 0;
        for (uint32_t p = c.begin; p < c.end; ++p) cap += db.off[L.occ[p] + 1] - db.off[L.occ[p]];
        child.reserve(c.end - c.begin, cap);
        uint32_t* out = child.items.data();
        size_t nt = 0;
        for (uint32_t p = c.begin; p < c.end; ++p) {
            const uint32_t t = L.occ[p];
            // items are sorted, so the ones after c.item start right after its position
            const uint32_t* b = db.items.data() + L.occpos[p] + 1;
            const uint32_t* e = db.items.data() + db.off[t + 1];
            uint32_t* start = out;
            uint64_t h = 0x9E3779B97F4A7C15ULL;
            for (; b != e; ++b) {
                const uint32_t j = *b;
                const bool keep = !((acc[j >> 6] >> (j & 63)) & 1ULL);   // drop patients now in the group
                *out = j;
                out += keep;
                if (keep) h = hash_step(h, j);
            }
            if (out == start) continue;           // nothing left to add from this atom
            std::memcpy(child.bits.data() + nt * W_, db.row(t), W_ * sizeof(uint64_t));
            child.weight[nt] = db.weight[t];
            child.hash[nt] = h;
            child.off[++nt] = uint32_t(out - child.items.data());
        }
        child.T = nt;
        child.NI = size_t(out - child.items.data());
        if (nt > 1) merge_duplicates(child, w.table);
        return child.T;
    }

    // Progress: each node's share of the whole search is split equally among
    // its candidate patients and credited when that candidate's sub-tree is
    // finished. Sub-trees differ enormously in size, so this is only a rough
    // guide: on test data it ran behind the true fraction by up to about 2x
    // (weighting by predicted sub-tree size was tried and did worse).
    void process(Worker& w, const uint64_t* G, int core, const DB& db, size_t depth, double share) {
        if (w.levels.size() <= depth) w.levels.resize(depth + 1);
        Worker::Level& L = w.levels[depth];
        prepare(w, db, L);
        const size_t nc = L.cands.size();
        if (nc == 0) { credit(w, share); return; }
        const double sub = share / double(nc);
        const int maxk = opt_.max_patients;
        for (size_t ci = 0; ci < nc; ++ci) {
            const Cand c = L.cands[ci];
            const int k = closure(G, db, L, c);
            if (k < 0) { credit(w, sub); continue; }
            emit(L.G.data(), c.sup, L.acc.data(), w.stats, w.buf);
            if ((maxk > 0 && k == maxk) || build_child(w, db, L, c) == 0) { credit(w, sub); continue; }
            const DB& child = L.child;
            if (idle_.load(std::memory_order_relaxed) > 0 && child.NI > 64) {
                // another thread is idle: give it this sub-tree
                auto task = std::make_unique<Task>();
                task->G = L.G;
                task->core = int(c.item);
                task->db = child.compact();
                task->share = sub;
                esync::lock_guard<esync::mutex> lock(qm_);
                ++pending_;
                queue_.push_back(std::move(task));
                qcv_.notify_one();
            } else {
                process(w, L.G.data(), int(c.item), child, depth + 1, sub);
            }
        }
        (void)core;
    }

public:
    // Knuth's (1975) estimator of the number of groups without enumerating them:
    // walk from the root choosing a random child at each level; the sum over
    // levels of the product of branching factors is an unbiased estimate of
    // the number of nodes (groups) in the search tree.
    std::vector<double> estimate(const std::vector<uint64_t>& G0, const DB& root, long walks, uint64_t seed) {
        Worker w;
        w.stats.init(m_);
        w.cnt.assign(m_, 0);
        w.wsum.assign(m_, 0);
        w.fill.assign(m_, 0);
        std::vector<double> est(walks);
        std::vector<size_t> ok, ok0;
        uint64_t state = seed * 0x9E3779B97F4A7C15ULL + 1;
        auto rnd = [&](size_t n) {   // splitmix64
            uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
            z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
            z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
            z ^= z >> 31;
            return size_t(z % n);
        };
        for (long walk = 0; walk < walks; ++walk) {
            const DB* db = &root;
            const uint64_t* G = G0.data();
            double prod = 1, sum = 0;
            for (size_t depth = 0;; ++depth) {
                if (w.levels.size() <= depth) w.levels.resize(depth + 1);
                Worker::Level& L = w.levels[depth];
                if (depth > 0 || walk == 0) {        // the root level is the same for every walk
                    prepare(w, *db, L);
                    ok.clear();
                    for (size_t ci = 0; ci < L.cands.size(); ++ci)
                        if (closure(G, *db, L, L.cands[ci]) >= 0) ok.push_back(ci);
                    if (depth == 0) ok0 = ok;
                } else {
                    ok = ok0;
                }
                if (ok.empty()) break;
                prod *= double(ok.size());
                sum += prod;
                const Cand c = L.cands[ok[rnd(ok.size())]];
                const int k = closure(G, *db, L, c);
                if (opt_.max_patients > 0 && k == opt_.max_patients) break;
                if (build_child(w, *db, L, c) == 0) break;
                db = &L.child;
                G = L.G.data();
            }
            est[walk] = sum;
            walks_done_.fetch_add(1, std::memory_order_relaxed);
        }
        return est;
    }

private:
    const Options& opt_;
    int m_, W_;
    std::vector<uint32_t> orig_;
    std::vector<uint64_t> use_mask_;
    GroupWriter* writer_;
    esync::mutex qm_;
    esync::condition_variable qcv_;
    std::deque<std::unique_ptr<Task>> queue_;
    std::atomic<int> idle_{0};
public:
    std::atomic<uint64_t> progress_{0};   // groups found so far (in steps of 4096)
    std::atomic<long> walks_done_{0};     // --estimate walks finished
private:
    std::vector<std::unique_ptr<std::atomic<double>>> done_;
public:
private:
    long pending_ = 0;
};

std::string fmt_count(uint64_t x) {
    std::string s = std::to_string(x);
    for (int p = int(s.size()) - 3; p > 0; p -= 3) s.insert(size_t(p), ",");
    return s;
}

std::string fmt_duration(double secs) {
    if (!(secs >= 0) || secs > 1e12) return "?";
    long t = long(secs + 0.5);
    std::ostringstream o;
    if (t >= 86400) o << t / 86400 << "d " << (t % 86400) / 3600 << "h";
    else if (t >= 3600) o << t / 3600 << "h " << (t % 3600) / 60 << "m";
    else if (t >= 60) o << t / 60 << "m " << t % 60 << "s";
    else o << t << "s";
    return o.str();
}

// Progress reporting: rewrites the --progress file every second (for
// elosa_run.m to display) and prints a line to stderr every 30 seconds.
// File format, one line: phase <tab> fraction done (-1 = unknown) <tab>
// groups so far <tab> seconds elapsed <tab> estimated seconds left (-1 = unknown)
class Monitor {
public:
    Monitor(const Options& o, std::chrono::steady_clock::time_point t0) : opt_(o), t0_(t0) {}
    ~Monitor() { stop(); }

    void phase(const std::string& name) { write(name, -1, 0, -1); }

    void start(const std::string& name, std::function<double()> fraction, std::function<uint64_t()> groups) {
        stop();
        done_ = false;
        const auto start = std::chrono::steady_clock::now();
        thread_ = std::make_unique<esync::thread>([this, name, fraction, groups, start] {
            auto last_msg = start;
            esync::unique_lock<esync::mutex> lock(m_);
            for (;;) {
                const bool finished = cv_.wait_for(lock, std::chrono::seconds(1), [&] { return done_; });
                if (finished) break;
                const auto now = std::chrono::steady_clock::now();
                const double f = fraction();
                const uint64_t g = groups();
                const double spent = std::chrono::duration<double>(now - start).count();
                const double left = (f > 1e-9 && spent >= 2) ? spent * (1 - f) / f : -1;
                write(name, f, g, left);
                if (!opt_.quiet && now - last_msg >= std::chrono::seconds(30)) {
                    last_msg = now;
                    std::cerr << "elosa: " << name << " " << std::fixed << std::setprecision(1) << 100 * f << "% done";
                    if (g) std::cerr << ", ~" << fmt_count(g) << " groups";
                    std::cerr << ", " << fmt_duration(elapsed()) << " so far";
                    if (left >= 0) std::cerr << ", about " << fmt_duration(left) << " left (rough)";
                    std::cerr << "\n" << std::defaultfloat;
                }
            }
        });
    }

    void stop() {
        if (!thread_) return;
        {
            esync::lock_guard<esync::mutex> lock(m_);
            done_ = true;
        }
        cv_.notify_all();
        thread_->join();
        thread_.reset();
    }

private:
    double elapsed() const {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0_).count();
    }
    void write(const std::string& name, double f, uint64_t g, double left) {
        if (opt_.progress.empty()) return;
        std::FILE* fp = std::fopen(opt_.progress.c_str(), "w");
        if (!fp) return;
        std::fprintf(fp, "%s\t%.6f\t%llu\t%.1f\t%.1f\n", name.c_str(), f, (unsigned long long)g, elapsed(), left);
        std::fclose(fp);
    }

    const Options& opt_;
    std::chrono::steady_clock::time_point t0_;
    esync::mutex m_;
    esync::condition_variable cv_;
    bool done_ = false;
    std::unique_ptr<esync::thread> thread_;
};

// Ends the program if the parent process (MATLAB) goes away, so a killed or
// crashed MATLAB session does not leave elosa running for hours.
class Watchdog {
public:
    explicit Watchdog(unsigned long pid) {
        if (!pid) return;
        thread_ = std::make_unique<esync::thread>([this, pid] {
            esync::unique_lock<esync::mutex> lock(m_);
            while (!cv_.wait_for(lock, std::chrono::seconds(1), [&] { return done_; })) {
                if (!process_alive(pid)) {
                    std::fprintf(stderr, "elosa: parent process %lu has ended; stopping\n", pid);
                    std::fflush(stderr);
                    std::_Exit(3);
                }
            }
        });
    }
    ~Watchdog() {
        if (!thread_) return;
        {
            esync::lock_guard<esync::mutex> lock(m_);
            done_ = true;
        }
        cv_.notify_all();
        thread_->join();
    }
private:
    esync::mutex m_;
    esync::condition_variable cv_;
    bool done_ = false;
    std::unique_ptr<esync::thread> thread_;
};

}  // namespace

int main(int argc, char** argv) {
    Options opt = parse_args(argc, argv);
    auto t0 = std::chrono::steady_clock::now();
    Watchdog watchdog(opt.parent_pid);
    Monitor monitor(opt, t0);
    try {
        monitor.phase("reading");
        Input in = read_input(opt.input);
        monitor.phase("preparing");
        const uint32_t n = in.n;

        // ---- patient order: grouped patients first (sorted by lesion size), then the rest
        std::vector<uint64_t> vox(n, 0);
        for (uint32_t v = 0; v < in.V; ++v) {
            const uint8_t* col = in.lesion.data() + size_t(v) * n;
            for (uint32_t p = 0; p < n; ++p) vox[p] += col[p] != 0;
        }
        std::vector<uint32_t> useidx, other;
        for (uint32_t p = 0; p < n; ++p) (in.use[p] ? useidx : other).push_back(p);
        if (opt.order != "input") {
            bool asc = opt.order == "ascending";
            std::stable_sort(useidx.begin(), useidx.end(), [&](uint32_t a, uint32_t b) {
                return asc ? vox[a] < vox[b] : vox[a] > vox[b]; });
        }
        const int m = int(useidx.size());
        if (m == 0) throw std::runtime_error("no patients are marked for grouping");
        std::vector<uint32_t> orig;             // new index -> original 1-based index
        std::vector<uint32_t> newidx(n);
        for (uint32_t p : useidx) { newidx[p] = uint32_t(orig.size()); orig.push_back(p + 1); }
        for (uint32_t p : other) { newidx[p] = uint32_t(orig.size()); orig.push_back(p + 1); }
        const int W = int((n + 63) / 64);

        // ---- atoms: distinct voxel patterns that involve at least one grouped patient
        std::vector<uint64_t> vbits;
        std::vector<uint32_t> vlist;
        vbits.reserve(size_t(in.V) * W / 4 + W);
        std::vector<uint64_t> row(W);
        for (uint32_t v = 0; v < in.V; ++v) {
            const uint8_t* col = in.lesion.data() + size_t(v) * n;
            std::fill(row.begin(), row.end(), 0);
            bool any_use = false;
            for (uint32_t p = 0; p < n; ++p)
                if (col[p]) {
                    uint32_t q = newidx[p];
                    row[q >> 6] |= 1ULL << (q & 63);
                    any_use |= q < uint32_t(m);
                }
            if (!any_use) continue;
            vlist.push_back(uint32_t(vbits.size() / W));
            vbits.insert(vbits.end(), row.begin(), row.end());
        }
        in.lesion.clear();
        in.lesion.shrink_to_fit();
        std::vector<uint32_t> ord(vlist.size());
        std::iota(ord.begin(), ord.end(), 0);
        std::sort(ord.begin(), ord.end(), [&](uint32_t a, uint32_t b) {
            return std::lexicographical_compare(vbits.begin() + size_t(a) * W, vbits.begin() + size_t(a + 1) * W,
                                                vbits.begin() + size_t(b) * W, vbits.begin() + size_t(b + 1) * W);
        });
        DB root;
        root.W = W;
        std::vector<uint64_t> rootacc(W, ~0ULL);
        uint64_t total_vox = 0;
        size_t natoms = 0;
        for (size_t k = 0; k < ord.size();) {
            size_t k2 = k + 1;
            const uint64_t* a = vbits.data() + size_t(ord[k]) * W;
            while (k2 < ord.size() && std::equal(a, a + W, vbits.data() + size_t(ord[k2]) * W)) ++k2;
            ++natoms;
            total_vox += k2 - k;
            for (int x = 0; x < W; ++x) rootacc[x] &= a[x];
            root.bits.insert(root.bits.end(), a, a + W);
            root.weight.push_back(k2 - k);
            for (int q = 0; q < m; ++q)
                if ((a[q >> 6] >> (q & 63)) & 1ULL) root.items.push_back(uint32_t(q));
            root.off.push_back(uint32_t(root.items.size()));
            k = k2;
        }
        vbits.clear();
        vbits.shrink_to_fit();
        root.T = root.weight.size();
        root.NI = root.items.size();

        // the root closure: grouped patients lesioned in every such voxel (usually none)
        std::vector<uint64_t> G0(W, 0);
        int k0 = 0;
        for (int q = 0; q < m; ++q)
            if (natoms && ((rootacc[q >> 6] >> (q & 63)) & 1ULL)) { G0[q >> 6] |= 1ULL << (q & 63); ++k0; }
        if (k0) {   // drop them from the item lists
            DB r2; r2.W = W;
            for (size_t t = 0; t < root.size(); ++t) {
                size_t before = r2.items.size();
                for (uint32_t p = root.off[t]; p < root.off[t + 1]; ++p) {
                    uint32_t q = root.items[p];
                    if (!((G0[q >> 6] >> (q & 63)) & 1ULL)) r2.items.push_back(q);
                }
                if (r2.items.size() == before) continue;
                r2.bits.insert(r2.bits.end(), root.row(t), root.row(t) + W);
                r2.weight.push_back(root.weight[t]);
                r2.off.push_back(uint32_t(r2.items.size()));
            }
            r2.T = r2.weight.size();
            r2.NI = r2.items.size();
            std::swap(root, r2);
        }

        int threads = opt.threads > 0 ? opt.threads : int(std::max(1u, esync::thread::hardware_concurrency()));
        if (!opt.quiet)
            std::cerr << "elosa: " << n << " patients (" << m << " grouped), " << total_vox << " lesioned voxels, "
                      << natoms << " distinct voxel patterns; " << threads << " threads\n";

        std::unique_ptr<GroupWriter> writer;
        if (!opt.out.empty()) writer = std::make_unique<GroupWriter>(opt.out, n, (opt.covering ? 1u : 0u) | (n < 65536 ? 2u : 0u));
        Search search(opt, m, W, orig, writer.get(), threads);
        Stats total;
        Stats rootstats; rootstats.init(m);
        OutBuf rootbuf;
        if (k0 && total_vox >= opt.min_voxels && (opt.max_patients == 0 || k0 <= opt.max_patients)) search.emit_root(G0, total_vox, rootacc.data(), rootstats, rootbuf);
        if (writer && !rootbuf.head.empty()) writer->write(rootbuf);
        if (opt.estimate > 0) {
            // estimates from each thread's walks, pooled
            std::vector<std::vector<double>> parts(threads);
            monitor.start("estimating", [&] { return double(search.walks_done_.load()) / double(opt.estimate); },
                          [] { return uint64_t(0); });
            std::vector<esync::thread> pool;
            for (int t = 0; t < threads; ++t) {
                long nw = opt.estimate / threads + (t < opt.estimate % threads ? 1 : 0);
                pool.emplace_back([&, t, nw] { parts[t] = search.estimate(G0, root, nw, opt.seed * 1000003ULL + t); });
            }
            for (auto& th : pool) th.join();
            monitor.stop();
            double sum = 0, sum2 = 0; long cnt = 0;
            for (auto& v : parts) for (double x : v) { sum += x; sum2 += x * x; ++cnt; }
            const double mean = sum / cnt + double(rootstats.groups);
            const double se = cnt > 1 ? std::sqrt(std::max(0.0, (sum2 - sum * sum / cnt) / (cnt - 1)) / cnt) : 0;
            double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            std::ostringstream s;
            s << "patients\t" << n << "\ngrouped_patients\t" << m << "\nlesioned_voxels\t" << total_vox
              << "\ndistinct_voxel_patterns\t" << natoms << "\nmin_voxels\t" << opt.min_voxels
              << "\nmax_patients\t" << opt.max_patients << "\nestimated_groups\t" << mean
              << "\nstandard_error\t" << se << "\nwalks\t" << cnt << "\nseconds\t" << secs << "\n";
            if (!opt.summary.empty()) {
                std::ofstream f(opt.summary);
                if (!f) throw std::runtime_error("cannot write " + opt.summary);
                f << s.str();
            }
            if (!opt.quiet) std::cerr << "elosa: estimated " << mean << " groups (standard error " << se << ", "
                                      << cnt << " walks) in " << secs << " s\n";
            std::cout << mean << " " << se << "\n";
            monitor.phase("done");
            return 0;
        }

        monitor.start("searching", [&] { return search.fraction(); }, [&] { return search.progress_.load(); });
        Task task{G0, -1, std::move(root), 1.0};
        search.run(std::move(task), threads, total);
        monitor.stop();
        monitor.phase("finishing");
        total.add(rootstats);
        if (writer && writer->failed()) throw std::runtime_error("writing " + opt.out + " failed (disk full?)");
        writer.reset();

        double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::ostringstream s;
        s << "patients\t" << n << "\ngrouped_patients\t" << m << "\nlesioned_voxels\t" << total_vox
          << "\ndistinct_voxel_patterns\t" << natoms << "\nmin_voxels\t" << opt.min_voxels
          << "\nmax_patients\t" << opt.max_patients << "\ngroups\t" << total.groups
          << "\nseconds\t" << secs << "\nthreads\t" << threads << "\n\n# groups by number of patients\npatients\tgroups\n";
        for (int k = 0; k <= m; ++k)
            if (total.by_size[k]) s << k << "\t" << total.by_size[k] << "\n";
        s << "\n# groups by region size (voxels in [from, to])\nfrom\tto\tgroups\n";
        for (int b = 0; b < 64; ++b)
            if (total.by_region_log2[b])
                s << (1ULL << b) << "\t" << ((b == 63) ? ~0ULL : (1ULL << (b + 1)) - 1) << "\t" << total.by_region_log2[b] << "\n";
        if (!opt.summary.empty()) {
            std::ofstream f(opt.summary);
            if (!f) throw std::runtime_error("cannot write " + opt.summary);
            f << s.str();
        }
        if (!opt.quiet) std::cerr << "elosa: " << fmt_count(total.groups) << " groups in " << secs << " s\n";
        std::cout << total.groups << "\n";
        monitor.phase("done");
    } catch (const std::exception& e) {
        std::cerr << "elosa: error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
