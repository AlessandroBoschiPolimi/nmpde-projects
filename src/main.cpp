#include <NeoHooke.hpp>

int main(int argc, char *argv[])
{
    using namespace pde;
    const unsigned int dim = NeoHooke::dim;
    const unsigned int r = 1;
    const unsigned int num_cells = 2;
    std::function<Point<dim> (const Point<dim> &)> h;
    NeoHooke problem = NeoHooke("", r, h, num_cells);
    problem.setup();
    problem.solve();
    return 0;
}
