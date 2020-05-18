// simple lock-free triple buffer
// unfair, but since operations should take a long time
// this shouldn't be a major issue.

#include "interface.h"

#include <atomic>

#include <cstdio>
#include <cstdlib>

#include <stdint.h>

#ifdef ASSERT
# undef ASSERT
#endif

#ifndef _NDEBUG

// provided for ease of debugging
void assertion_failure() {
    std::fprintf(stderr, "ASSERTION FAILURE\n");
}

static void sync_assert(bool cond, const char *message) {
    if (!cond) {
        assertion_failure();
        std::fprintf(stderr, message);
        std::abort();
    }
}

# define ASSERT(c) sync_assert(c, #c)
#else
# define ASSERT(c)
#endif

union OwnedBufferValue {
    uint64_t value64;
    struct {
        int32_t value32[NUM_OWNERS];
    } as_arr;
};

TripleBuffer<Grid> grids;
TripleBuffer<SimParams> param_buf;

template<typename T>
T *TripleBuffer<T>::get_init() {
    return &buffers[0];
}

template<typename T>
void TripleBuffer<T>::init() {
    for (auto &out : buffers) {
        out = buffers[0];
    }
}

inline static int32_t other_buf(OwnedBufferValue val) {
    constexpr int SUM_INDICES = 0 + 1 + 2;
    ASSERT(val.as_arr.value32[READER] >= 0 &&
           val.as_arr.value32[READER] <= 2);
    ASSERT(val.as_arr.value32[WRITER] >= 0 &&
           val.as_arr.value32[WRITER] <= 2);
    return SUM_INDICES
        - val.as_arr.value32[READER]
        - val.as_arr.value32[WRITER];
}

template<typename T>
T *TripleBuffer<T>::get_current(TripleBuffer<T>::Owner owner) {
    ASSERT(owner == READER || owner == WRITER);
    OwnedBufferValue val;
    uint64_t read;
    read = val.value64 = owned_buffers.load();
    return &buffers[val.as_arr.value32[WRITER]];
}

template<typename T>
T *TripleBuffer<T>::swap(TripleBuffer<T>::Owner owner) {
    ASSERT(owner == READER || owner == WRITER);
    OwnedBufferValue val;
    uint64_t read;
    int32_t ret;
    while (true) {
        read = val.value64 = owned_buffers.load();
        ret = val.as_arr.value32[owner] = other_buf(val);
        if (owned_buffers.compare_exchange_weak(read, val.value64)) {
            break;
        }
    }
    ASSERT(ret >= 0 && ret <= 2);
    return &buffers[ret];
}

template class TripleBuffer<Grid>;
template class TripleBuffer<SimParams>;
