// simple lock-free triple buffer
// unfair, but since operations should take a long time
// this shouldn't be a major issue.

#include "interface.h"

#include <atomic>

#include <cstdio>
#include <cstdlib>

#include <stdint.h>

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

# ifdef ASSERT
#  undef ASSERT
# endif
# define ASSERT(c) sync_assert(c, #c)
#else
# define ASSERT(c)
#endif

constexpr int NUM_BUFFERS = 3;

union OwnedGridValue {
    uint64_t value64;
    struct {
        int32_t value32[NUM_OWNERS];
    } as_arr;
};

static std::atomic_uint64_t owned_grids(0x000000000000001UL);

static Grid grids[NUM_BUFFERS];

Grid *get_init_grid() {
    return &grids[0];
}

void init_grids() {
    grids[0].updated = false;
    for (auto &out : grids) {
        out = grids[0];
    }
}

inline static int32_t other_grid(OwnedGridValue val) {
    constexpr int SUM_GRID_INDICES = 0 + 1 + 2;
    ASSERT(val.as_arr.value32[READER] >= 0 &&
           val.as_arr.value32[READER] <= 2);
    ASSERT(val.as_arr.value32[WRITER] >= 0 &&
           val.as_arr.value32[WRITER] <= 2);
    return SUM_GRID_INDICES
        - val.as_arr.value32[READER]
        - val.as_arr.value32[WRITER];
}

Grid *get_current_grid(GridOwner owner) {
    ASSERT(owner == READER || owner == WRITER);
    OwnedGridValue val;
    uint64_t read;
    read = val.value64 = owned_grids.load();
    return &grids[val.as_arr.value32[WRITER]];
}

Grid *grid_swap(GridOwner owner) {
    ASSERT(owner == READER || owner == WRITER);
    // if writer, mark as updated so reader knows
    // otherwise, unmark before releasing
    OwnedGridValue val;
    uint64_t read;
    read = val.value64 = owned_grids.load();
    grids[val.as_arr.value32[WRITER]].updated =
        owner == WRITER;
    int32_t ret;
    while (true) {
        ret = val.as_arr.value32[owner] = other_grid(val);
        if (owned_grids.compare_exchange_weak(read, val.value64)) {
            break;
        }
        read = val.value64 = owned_grids.load();
    }
    return &grids[ret];
}
