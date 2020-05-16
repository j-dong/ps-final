#include "interface.h"

static void step(Grid *grid, const Grid *prev);

void sim_init_grid(Grid *grid) {
    // ...
    (void) grid;
}

void sim_main() {
    Grid *prev = get_current_grid(WRITER);
    while (running.load(std::memory_order_relaxed)) {
        Grid *next = grid_swap(WRITER);
        step(next, prev);
        prev = next;
    }
}

void step(Grid *grid, const Grid *prev) {
    // ...
    *grid = *prev;
}
