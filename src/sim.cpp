#include "interface.h"

#include <cmath>
#include <cstdio>

#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

using namespace Eigen;

static SimParams params;

static void step(Grid *grid, const Grid *prev);

void sim_init_grid(Grid *grid) {
    for (int y = 0; y < 20; y++) for (int x = 0; x < 20; x++) {
        grid->density[y + HEIGHT - 40][WIDTH / 2 - 10 + x] = 1.0;
        grid->temperature[y + HEIGHT - 40][WIDTH / 2 - 10 + x] = 20;
    }
}

void sim_main() {
    Grid *prev = grids.get_current(WRITER);
    while (running.load(std::memory_order_relaxed)) {
        Grid *next = grids.swap(WRITER);
        auto new_params = param_buf.swap(READER);
        if (new_params->updated) {
            params = *new_params;
            new_params->updated = false;
        }
        step(next, prev);
        next->updated = true;
        /*
        for (int y = 0; y < HEIGHT; y++) {
            for (int x = 0; x < WIDTH; x++) {
                std::cout << (next->temperature[y][x] > 10.0f ? 'A' : ' ');
            }
            std::cout << "\n";
        }
        std::cout << std::endl;
        */
        prev = next;
    }
}

template<int W, int H>
double interpolate(const double (&values)[H][W], Vector2d position) {
    if (position(0) <= 0.0) position(0) = 0.0;
    if (position(1) <= 0.0) position(1) = 0.0;
    if (position(0) > W - 1) position(0) = W - 1;
    if (position(1) > H - 1) position(1) = H - 1;
    int ix = (int) position(0),
        iy = (int) position(1);
    if (ix >= W - 1) ix = W - 2;
    if (iy >= H - 1) iy = H - 2;
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

void process_forces(Grid *grid) {
    for (int y = 0; y < HEIGHT; y++) for (int x = 0; x < WIDTH; x++) {
        Vector2d pt(x, y - 0.5);
        double density = interpolate(grid->density, pt);
        double temp = interpolate(grid->temperature, pt);
        grid->velocity_y[y][x] += params.timestep * (-params.alpha * density + params.beta * temp);
    }
    for (int y = 0; y < HEIGHT; y++) for (int x = 0; x < WIDTH; x++)
    {
      Vector2d vel_u = interpolateVelocity(*grid, Vector2d(x,y+1));
      Vector2d vel_d = interpolateVelocity(*grid, Vector2d(x,y-1));
      Vector2d vel_l = interpolateVelocity(*grid, Vector2d(x+1,y));
      Vector2d vel_r = interpolateVelocity(*grid, Vector2d(x-1,y));

      //Since 2-D we only need z direction of omega
      grid->vorticity[y][x] = 0.5 * (vel_r(1) - vel_l(1) - vel_u(0) + vel_d(1));
    }

    for (int y = 1; y < HEIGHT - 1; y++) for (int x = 1; x < WIDTH - 1; x++)
    {
      Vector2d Nx;
      //Used for normalization
      using std::abs;
      Nx(0) = abs(grid->vorticity[y][x]) - abs(grid->vorticity[y][x-1]);
      Nx(1) = 0.25 * (abs(grid->vorticity[y+1][x]) + abs(grid->vorticity[y+1][x-1]) - abs(grid->vorticity[y-1][x]) - abs(grid->vorticity[y-1][x-1]));
      Nx /= Nx.norm() + 1e-5;
      Vector2d Ny;
      //Used for normalization
      Ny(1) = abs(grid->vorticity[y][x]) - abs(grid->vorticity[y-1][x]);
      Ny(0) = 0.25 * (abs(grid->vorticity[y][x+1]) + (grid->vorticity[y-1][x+1]) - abs(grid->vorticity[y][x-1]) - abs(grid->vorticity[y-1][x-1]));
      Ny /= Ny.norm() + 1e-5;
      double omega_x = 0.5* (grid->vorticity[y][x] + grid->vorticity[y][x-1]);
      double omega_y = 0.5* (grid->vorticity[y][x] + grid->vorticity[y-1][x]);
      grid->velocity_x[y][x] += Nx(1) * omega_x * params.epsilon * params.timestep;
      grid->velocity_y[y][x] -= Ny(0) * omega_y * params.epsilon * params.timestep;
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
        b(y) /= params.timestep;
    }
    A.setFromTriplets(rows.begin(), rows.end());
    rows.clear();
    ConjugateGradient<Mat, Lower|Upper, IncompleteCholesky<double>> solver;
    // solver.setMaxIterations(20);
    // solver.setTolerance(1e-8);
    solver.compute(A);
    Map<VectorXd> pressureMap((double *) grid->pressure, N);
    VectorXd temp = solver.solve(b);
    pressureMap = temp;
}

void step(Grid *grid, const Grid *prev) {
    *grid = *prev;
    for (int y = 0; y < HEIGHT + 1; y++) for (int x = 0; x < WIDTH + 1; x++) {
        // update velocity to obtain u*
        grid->velocity_x[y][x] = interpolate(
            prev->velocity_x,
            Vector2d(x, y) - params.timestep *
                interpolateVelocity(*prev, Vector2d(x - 0.5, y))
        );
        grid->velocity_y[y][x] = interpolate(
            prev->velocity_y,
            Vector2d(x, y) - params.timestep *
                interpolateVelocity(*prev, Vector2d(x, y - 0.5))
        );
    }
    process_forces(grid);
    calculate_pressure(grid);
    for (int y = 1; y < HEIGHT; y++) for (int x = 0; x < WIDTH; x++) {
        // update velocity based on pressure
        grid->velocity_x[y][x] -= params.timestep * (grid->pressure[y][x] - grid->pressure[y][x - 1]);
        grid->velocity_y[y][x] -= params.timestep * (grid->pressure[y][x] - grid->pressure[y - 1][x]);
    }
    // advect temperature and density
    for (int y = 0; y < HEIGHT; y++) for (int x = 0; x < WIDTH; x++) {
        Vector2d pt = Vector2d(x, y) - params.timestep *
            interpolateVelocity(*grid, Vector2d(x, y));
        grid->density[y][x] = interpolate(prev->density, pt);
        grid->temperature[y][x] = interpolate(prev->temperature, pt);
    }
    for (int z = 0; z < 10; z++) {
        int y = 2;
        int x = (WIDTH - 10) / 2 + z;
        double d = grid->density[y][x];
        grid->density[y][x] += 1.0;
        grid->temperature[y][x] = (grid->temperature[y][x] * d + 1.0 * 100.0) / grid->density[y][x];
        // grid->velocity_y[y][x] = 1.0;
    }
}
