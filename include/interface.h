#ifndef __INCLUDED_INTERFACE_H__
#define __INCLUDED_INTERFACE_H__

#include <atomic>

extern std::atomic_bool running;

constexpr int WIDTH = 64, HEIGHT = 64;

enum BufferOwner {
    READER,
    WRITER,
    NUM_OWNERS,
};

template<typename T>
class TripleBuffer {
    T buffers[3];
    std::atomic_uint64_t owned_buffers = {0x00000000'00000001UL};

public:
    using Owner = BufferOwner;

    T *swap(Owner owner);
    T *get_current(Owner owner);

    T *get_init();
    void init();
};

struct alignas(4 * sizeof(double)) Grid
{
    double pressure[HEIGHT][WIDTH];
    double density[HEIGHT][WIDTH];
    double temperature[HEIGHT][WIDTH];
    double vorticity[HEIGHT][WIDTH];
    // velocity_x[y][x] connects (x - 1, y) with (x, y)
    double velocity_x[HEIGHT + 4][WIDTH + 4];
    // velocity_y[y][x] connects (x, y - 1) with (x, y)
    double velocity_y[HEIGHT + 4][WIDTH + 4];
    // feel free to add
    bool updated; // handled automatically for you
};

struct SimParams {
    double timestep = 0.01;
    double alpha = 0.2, beta = 1.0;
    double epsilon = 0.01;
    bool updated = false;
};

extern TripleBuffer<Grid> grids;
extern TripleBuffer<SimParams> param_buf;

extern template class TripleBuffer<Grid>;
extern template class TripleBuffer<SimParams>;

#endif
