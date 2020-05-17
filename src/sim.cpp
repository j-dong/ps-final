#include "interface.h"

#include <cmath>
#include <cstdio>

#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

using namespace Eigen;

static void step(Grid *grid, const Grid *prev);

void sim_init_grid(Grid *grid) {
    for (int y = 0; y < 20; y++) for (int x = 0; x < 20; x++) {
        grid->density[y + HEIGHT - 40][WIDTH / 2 - 10 + x] = 1.0;
        grid->temperature[y + HEIGHT - 40][WIDTH / 2 - 10 + x] = 20;
    }
    for (int y = 0; y < HEIGHT; y++) for (int x = 0; x < WIDTH; x++) {
        grid->velocity_y[y][x] = 1.0;
    }
}

void sim_main() {
    Grid *prev = get_current_grid(WRITER);
    while (running.load(std::memory_order_relaxed)) {
        Grid *next = grid_swap(WRITER);
        step(next, prev);
        next->updated = true;
        printf("updated %p\n", next);
        /*
        for (int y = 0; y < HEIGHT; y++) {
            for (int x = 0; x < WIDTH; x++) {
                std::cout << (next->temperature[y][x] > 10.0f ? 'A' : ' ');
            }
            std::cout << "\n";
        }
        */
        std::cout << std::endl;
        prev = next;
    }
}

template<int W, int H>
double interpolate(const double (&values)[H][W], Vector2d position) {
    if (position(0) < 0.0) position(0) = 0.0;
    if (position(1) < 0.0) position(1) = 0.0;
    if (position(0) > W) position(0) = W;
    if (position(1) > H) position(1) = H;
    int ix = (int) std::floor(position(0)),
        iy = (int) std::floor(position(1));
    if (ix >= W) ix = W - 1;
    if (iy >= H) iy = H - 1;
    double fx = position(0) - ix, fy = position(1) - iy;
    double vx0 = values[iy]    [ix] * (1.0 - fx) + values[iy]    [ix + 1] * fx;
    double vx1 = values[iy + 1][ix] * (1.0 - fx) + values[iy + 1][ix + 1] * fx;
    return vx0 * (1.0 - fy) + vx1 * fy;
}

// estimate velocity at position by interpolating nearest neighbors
Vector2d interpolateVelocity(const Grid &grid, Vector2d position) {
    double vx = interpolate(grid.velocity_x, position + Vector2d(0.5, 0.0));
    double vy = interpolate(grid.velocity_y, position + Vector2d(0.0, 0.5));
    return Vector2d(vx, vy);
}

constexpr double TIMESTEP = 0.01;
constexpr double ALPHA    = 0.2;
constexpr double BETA     = 1.0;

void process_forces(Grid *grid) {
    for (int y = 0; y < HEIGHT; y++) for (int x = 0; x < WIDTH; x++) {
        grid->velocity_y[y][x] += TIMESTEP * (-ALPHA * grid->density[y][x] + BETA * grid->temperature[y][x]);
    }
}

inline int INDEX(int x, int y) { return x + y * WIDTH; }

bool is_valid(int x, int y) {
    if (x < 0 || x >= WIDTH) return false;
    if (y < 0 || y >= HEIGHT) return false;
    // TODO: solid cells
    return true;
}

void calculate_pressure(Grid *grid) {
    typedef SparseMatrix<double> Mat;
    typedef VectorXd Vec;
    int N = WIDTH * HEIGHT;
    Mat A(N, N);
    Vec b(N);
    // fill in Laplacian matrix
    std::vector<Eigen::Triplet<double>> rows;
    rows.reserve(5 * N);
    int NEIGHBOR_OFFSETS[][2] = {
        {-1, 0},
        { 1, 0},
        { 0,-1},
        { 0, 1},
    };
    for (int y = 0; y < HEIGHT; y++) for (int x = 0; x < WIDTH; x++) {
        rows.emplace_back(INDEX(x, y), INDEX(x, y), -4.0);
        for (const int (&d)[2] : NEIGHBOR_OFFSETS) {
            int dx = x + d[0], dy = y + d[1];
            if (is_valid(dx, dy)) {
                rows.emplace_back(INDEX(x, y), INDEX(dx, dy), 1.0);
            }
        }
        b(y) = -grid->velocity_x[y][x] + grid->velocity_x[y][x + 1]
             + -grid->velocity_y[y][x] + grid->velocity_y[y + 1][x];
        b(y) /= TIMESTEP;
    }
    A.setFromTriplets(rows.begin(), rows.end());
    rows.clear();
    ConjugateGradient<Mat, Lower|Upper, IncompleteCholesky<double>> solver;
    solver.compute(A);
    Map<VectorXd> pressureMap((double *) grid->pressure, N);
    pressureMap = solver.solve(b);
}

void step(Grid *grid, const Grid *prev) {
    std::cout << "Processing next step" << std::endl;
    *grid = *prev;
    std::cout << "updating velocity due to advection" << std::endl;
    for (int y = 0; y < HEIGHT + 1; y++) for (int x = 0; x < WIDTH + 1; x++) {
        // update velocity to obtain u*
        grid->velocity_x[y][x] = interpolate(
            prev->velocity_x,
            Vector2d(x, y) - TIMESTEP *
                interpolateVelocity(*prev, Vector2d(x - 0.5, y))
        );
        grid->velocity_y[y][x] = interpolate(
            prev->velocity_y,
            Vector2d(x, y) - TIMESTEP *
                interpolateVelocity(*prev, Vector2d(x, y - 0.5))
        );
    }
    std::cout << "processing forces" << std::endl;
    process_forces(grid);
    std::cout << "calculating pressure" << std::endl;
    calculate_pressure(grid);
    std::cout << "updating velocity due to pressure" << std::endl;
    for (int y = 1; y < HEIGHT; y++) for (int x = 0; x < WIDTH; x++) {
        // update velocity based on pressure
        grid->velocity_x[y][x] -= TIMESTEP * (grid->pressure[y][x] - grid->pressure[y][x - 1]);
        grid->velocity_y[y][x] -= TIMESTEP * (grid->pressure[y][x] - grid->pressure[y - 1][x]);
    }
    std::cout << "updating temperature and density" << std::endl;
    // advect temperature and density
    for (int y = 0; y < HEIGHT; y++) for (int x = 0; x < WIDTH; x++) {
        Vector2d pt = Vector2d(x, y) - TIMESTEP *
            interpolateVelocity(*grid, Vector2d(x, y));
        grid->density[y][x] = interpolate(prev->density, pt);
        grid->temperature[y][x] = interpolate(prev->temperature, pt);
    }
    for (int z = 0; z < 10; z++) {
        int y = 1;
        int x = (WIDTH - 10) / 2 + z;
        double d = grid->density[y][x];
        grid->density[y][x] += 1.0;
        grid->temperature[y][x] = (grid->temperature[y][x] * d + 1.0 * 100.0) / grid->density[y][x];
        grid->velocity_y[y][x] = 1.0;
    }
}
