#ifndef __INCLUDED_INTERFACE_H__
#define __INCLUDED_INTERFACE_H__

#include <atomic>

extern std::atomic_bool running;

constexpr int WIDTH = 300, HEIGHT = 200;

enum GridOwner {
    READER, // read-only
    WRITER, // read/write
    NUM_OWNERS,
};

struct alignas(4 * sizeof(double)) Grid {
    double pressure[HEIGHT][WIDTH];
    double density[HEIGHT][WIDTH];
    double temperature[HEIGHT][WIDTH];
    // velocity_horizontal[y][x]
    // connects (x - 1, y) with (x, y)
    double velocity_horizontal[HEIGHT + 4][WIDTH + 4][2];
    // velocity_vertical[y][x]
    // connects (x, y - 1) with (x, y)
    double velocity_vertical[HEIGHT + 4][WIDTH + 4][2];
    // feel free to add
    bool updated; // handled automatically for you
};

// note: for writer, safe to read from previous pointer
Grid *grid_swap(GridOwner owner);

Grid *get_current_grid(GridOwner owner);

Grid *get_init_grid();
// initialize all grids with given data
void init_grids();

#endif
