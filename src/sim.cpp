#include "interface.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include <fstream>
#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>

constexpr bool DEBUG = false;

using namespace Eigen;

constexpr int N = WIDTH * HEIGHT;

void print_total_velocities(int line, Grid *grid) {
    double total_sum_x = 0.0;
    double total_sum_y = 0.0;
    for (int y = 0; y < HEIGHT + 1; y++) for (int x = 0; x < WIDTH + 1; x++) {
        total_sum_x += std::abs(grid->velocity_x[y][x]);
        total_sum_y += std::abs(grid->velocity_y[y][x]);
    }
    if (DEBUG) std::cout << "total at line " << line << ": " << total_sum_x << ", " << total_sum_y << std::endl;
}

#define PV do { if (DEBUG) print_total_velocities(__LINE__, grid); } while (false)

static SimParams params;
static SparseMatrix<double> laplacian;
static ConjugateGradient<SparseMatrix<double>, Lower|Upper, IncompleteCholesky<double>> solver;
// static SimplicialLDLT<SparseMatrix<double>> solver;

static void step(Grid *grid, const Grid *prev);

void sim_init_grid(Grid *grid) {
    std::memset(grid, 0, sizeof *grid);
    for (int y = 0; y < HEIGHT + 1; y++) for (int x = 0; x < WIDTH + 1; x++) {
        grid->velocity_x[y][x] = 0.0;
        grid->velocity_y[y][x] = 0.0;
    }
    for (int y = 0; y < 20; y++) for (int x = 0; x < 20; x++) {
        grid->density[y + HEIGHT - 40][WIDTH / 2 - 10 + x] = 1.0;
        grid->temperature[y + HEIGHT - 40][WIDTH / 2 - 10 + x] = 20;
    }
}

inline int INDEX(int x, int y) { return x + y * WIDTH; }

bool is_valid(int x, int y) {
    if (x < 0 || x >= WIDTH) return false;
    if (y < 0 || y >= HEIGHT) return false;
    // TODO: solid cells
    return true;
}

void sim_init() {
    solver.setMaxIterations(40);
    solver.setTolerance(1e-10);
    std::vector<Eigen::Triplet<double>> rows;
    rows.reserve(5 * N);
    // fill in Laplacian matrix
    int NEIGHBOR_OFFSETS[][2] = {
        {-1, 0},
        { 1, 0},
        { 0,-1},
        { 0, 1},
    };
    for (int y = 0; y < HEIGHT; y++) for (int x = 0; x < WIDTH; x++) {
        int neighbor_count = 0;
        for (const int (&d)[2] : NEIGHBOR_OFFSETS) {
            int dx = x + d[0], dy = y + d[1];
            if (is_valid(dx, dy)) {
                rows.emplace_back(INDEX(x, y), INDEX(dx, dy), -1.0);
                neighbor_count += 1;
            }
        }
        // rows.emplace_back(INDEX(x, y), INDEX(x, y), neighbor_count);
        rows.emplace_back(INDEX(x, y), INDEX(x, y), neighbor_count);
    }
    laplacian.resize(N, N);
    laplacian.setFromTriplets(rows.begin(), rows.end());
    rows.clear();
    solver.compute(laplacian);
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
double querySafe(const double (&values)[H][W], int x, int y) {
    if (x < 0 || y < 0 || x >= W || y >= H) return 0.0;
    return values[y][x];
}

template<int W, int H>
double interpolate(const double (&values)[H][W], Vector2d position) {
    int ix = (int) std::floor(position(0)),
        iy = (int) std::floor(position(1));
    double fx = position(0) - ix, fy = position(1) - iy;
    double vx0 = querySafe(values, ix, iy  ) * (1.0 - fx) + querySafe(values, ix+1, iy  ) * fx;
    double vx1 = querySafe(values, ix, iy+1) * (1.0 - fx) + querySafe(values, ix+1, iy+1) * fx;
    return vx0 * (1.0 - fy) + vx1 * fy;
}

// estimate velocity at position by interpolating nearest neighbors
Vector2d interpolateVelocity(const Grid &grid, Vector2d position) {
    double vx = interpolate(grid.velocity_x, position + Vector2d(0.5, 0.0));
    double vy = interpolate(grid.velocity_y, position + Vector2d(0.0, 0.5));
    return Vector2d(vx, vy);
}

void process_forces(Grid *grid) {
    memset(grid->force_x, 0, sizeof grid->force_x);
    memset(grid->force_y, 0, sizeof grid->force_y);
    for (int y = 0; y < HEIGHT; y++) for (int x = 0; x < WIDTH; x++) {
        grid->force_y[y][x] += -params.alpha * grid->density[y][x] + params.beta * grid->temperature[y][x];
    }
    for (int y = 0; y < HEIGHT; y++) for (int x = 0; x < WIDTH; x++)
    {
      Vector2d vel_u = interpolateVelocity(*grid, Vector2d(x,y+1));
      Vector2d vel_d = interpolateVelocity(*grid, Vector2d(x,y-1));
      Vector2d vel_l = interpolateVelocity(*grid, Vector2d(x+1,y));
      Vector2d vel_r = interpolateVelocity(*grid, Vector2d(x-1,y));

      //Since 2-D we only need z direction of omega
      grid->vorticity[y][x] = 0.5 * (vel_r(1) - vel_l(1) - vel_u(0) + vel_d(0));
    }
    for (int y = 1; y < HEIGHT - 1; y++) for (int x = 1; x < WIDTH - 1; x++)
    {
      Vector2d N;
      //Used for normalization
      using std::abs;
      N(0) = 0.5 * (abs(grid->vorticity[y][x+1]) - abs(grid->vorticity[y][x-1]));
      N(1) = 0.5 * (abs(grid->vorticity[y+1][x]) - abs(grid->vorticity[y-1][x]));
      N /= N.norm() + 1e-5;
      grid->force_x[y][x] += N(1) * grid->vorticity[y][x] * params.epsilon;
      grid->force_y[y][x] -= N(0) * grid->vorticity[y][x] * params.epsilon;
    }
}

void apply_force(Grid *grid) {
    for (int y = 0; y < HEIGHT + 1; y++) for (int x = 0; x < WIDTH + 1; x++) {
        grid->velocity_x[y][x] += params.timestep * interpolate(grid->force_x, Vector2d(x - 0.5, y));
        grid->velocity_y[y][x] += params.timestep * interpolate(grid->force_y, Vector2d(x, y - 0.5));
    }
}

void calculate_pressure(Grid *grid) {
    VectorXd b(N);
    double sum_b = 0.0;
    double dot_b = 0.0;
    for (int y = 0; y < HEIGHT; y++) for (int x = 0; x < WIDTH; x++) {
        int n = INDEX(x, y);
        b(n) = -grid->velocity_x[y][x] + grid->velocity_x[y][x + 1]
             + -grid->velocity_y[y][x] + grid->velocity_y[y + 1][x];
        b(n) *= -1.0;
        sum_b += std::abs(b(n));
        dot_b += b(n);
    }
    if (DEBUG) std::cout << "sum of b's: " << sum_b << std::endl;
    if (DEBUG) std::cout << "dot b with 1: " << dot_b << std::endl;
    Map<VectorXd> pressureMap((double *) grid->pressure, N);
    VectorXd temp = solver.solve(b);
    if (DEBUG) std::cout << "solver error: " << ((laplacian * temp) - b).norm() << std::endl;
    // if (b.norm() < 1e-10) {
    //     temp.setZero();
    // }
    double sum_p = 0.0;
    for (int i = 0; i < N; i++) {
        sum_p += std::abs(temp(i));
    }
    if (DEBUG) std::cout << "sum of P: " << sum_p << std::endl;
    pressureMap = temp;
}

bool have_exported = false;

void step(Grid *grid, const Grid *prev) {
    *grid = *prev;
    PV;
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
    PV;
    apply_force(grid);
    PV;
    calculate_pressure(grid);
    {
        double sum_div = 0.0;
        for (int y = 1; y < HEIGHT - 1; y++) for (int x = 1; x < WIDTH - 1; x++) {
            double div = -grid->velocity_x[y][x] + grid->velocity_x[y][x+1]
                         -grid->velocity_y[y][x] + grid->velocity_y[y+1][x];
            sum_div += std::abs(div);
            // if (std::fabs(div) > 1e-2) std::cout << "+++" << div << std::endl;
        }
        if (DEBUG) std::cout << "sum of div before: " << sum_div << std::endl;
    }
    for (int y = 1; y < HEIGHT; y++) for (int x = 1; x < WIDTH; x++) {
        // update velocity based on pressure
        grid->velocity_x[y][x] -= (grid->pressure[y][x] - grid->pressure[y][x - 1]);
        grid->velocity_y[y][x] -= (grid->pressure[y][x] - grid->pressure[y - 1][x]);
    }
    PV;
    {
        double sum_div = 0.0;
        for (int y = 1; y < HEIGHT - 1; y++) for (int x = 1; x < WIDTH - 1; x++) {
            double div = -grid->velocity_x[y][x] + grid->velocity_x[y][x+1]
                         -grid->velocity_y[y][x] + grid->velocity_y[y+1][x];
            sum_div += std::abs(div);
            // if (std::fabs(div) > 1e-2) std::cout << "+++" << div << std::endl;
        }
        if (DEBUG) std::cout << "sum of div after: " << sum_div << std::endl;
    }
    for (int i = 0; i <= WIDTH; i++) {
        grid->velocity_y[0][i] = 0;
        grid->velocity_y[HEIGHT][i] = 0;
    }
    for (int i = 0; i <= HEIGHT; i++) {
        grid->velocity_y[i][0] = 0;
        grid->velocity_y[i][WIDTH] = 0;
    }
    // advect temperature and density
    for (int y = 0; y < HEIGHT; y++) for (int x = 0; x < WIDTH; x++) {
        Vector2d pt = Vector2d(x, y) - params.timestep *
            interpolateVelocity(*grid, Vector2d(x, y));
        grid->density[y][x] = interpolate(prev->density, pt);
        if (grid->density[y][x] < 0.0) {
            std::cout << "NEGATIVE DENSITY!" << std::endl;
        }
        grid->temperature[y][x] = interpolate(prev->temperature, pt);
        // heat transfer
        double old = grid->temperature[y][x];
        grid->temperature[y][x] -= params.kappa * grid->temperature[y][x] * params.timestep;
        if (old * grid->temperature[y][x] < 0.0) {
            grid->temperature[y][x] = 0.0;
        }
    }
    PV;
    if (params.emitter_density > 1e-6) {
        for (int z = 0; z < 10; z++) {
            int y = 2;
            int x = (WIDTH - 10) / 2 + z;
            double d = grid->density[y][x];
            double dens = params.emitter_density;
            grid->density[y][x] += dens;
            grid->temperature[y][x] = (grid->temperature[y][x] * d + dens * params.emitter_temp) / grid->density[y][x];
            // grid->velocity_y[y][x] = 1.0;
        }
    }
    PV;
    if (params.want_to_export && !have_exported) {
        have_exported = true;
        std::ofstream f("out.csv");
        for (int y = 0; y <= HEIGHT; y++) for (int x = 0; x <= WIDTH; x++) {
            f << x << "," << y << "," << grid->velocity_x[y][x] << "," << grid->velocity_y[y][x] << std::endl;
        }
    }
}
