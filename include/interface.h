#ifndef __INCLUDED_INTERFACE_H__
#define __INCLUDED_INTERFACE_H__

#include <atomic>

extern std::atomic_bool running;

constexpr int WIDTH = 128, HEIGHT = 128;

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
    double force_x[HEIGHT][WIDTH];
    double force_y[HEIGHT][WIDTH];
    // velocity_x[y][x] connects (x - 1, y) with (x, y)
    double velocity_x[HEIGHT + 4][WIDTH + 4];
    // velocity_y[y][x] connects (x, y - 1) with (x, y)
    double velocity_y[HEIGHT + 4][WIDTH + 4];
    // feel free to add
    bool updated; // handled automatically for you
};

struct SimParams {
    double timestep = 0.01;
    // double alpha = 0.0, beta = 0.0;
    double alpha = 2.0, beta = 10.0;
    double epsilon = 20.0;
    double emitter_density = 0.1;
    double emitter_temp = 50.0;
    bool want_to_export = false;
    bool updated = false;
};

extern TripleBuffer<Grid> grids;
extern TripleBuffer<SimParams> param_buf;

extern template class TripleBuffer<Grid>;
extern template class TripleBuffer<SimParams>;

#endif
