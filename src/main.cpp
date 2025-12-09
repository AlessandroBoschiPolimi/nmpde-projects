#include <NeoHooke.hpp>

int main(/*int argc, char *argv[]*/)
{
    using namespace pde;
    const unsigned int dim = MechanicalDisplacement::dim;
    const unsigned int r = 1;

    //Material constant
    const double C = 3; //Pa
    //Number of cells inside a hypercube
    const unsigned int num_cells = 8;
    const double tau_0 = 0.5;
    //Pulling the cube appart <-|  |-> in x,y directions
    //The dirichlet conditions are set to 0 on the z boundaries
    //The object can change volume, since lambda = 0 in our model
    const auto h  = [&tau_0](const Point<dim> &p) {
	double small_tol = 1e-13;
	if (std::abs(p[0]) < small_tol && std::abs(p[1]) > small_tol)
		return Point<dim>(-tau_0, 0, 0);
	else if (std::abs(p[0]) > (1 - small_tol) && std::abs(p[1]) > small_tol)
		return Point<dim>(tau_0, 0, 0);
	else if (std::abs(p[1]) < small_tol && std::abs(p[0]) > small_tol)
		return Point<dim>(0, -tau_0, 0);
	else if (std::abs(p[1]) > (1 - small_tol) && std::abs(p[0]) > small_tol)
		return Point<dim>(0, tau_0, 0);
	else
		return Point<dim>(0, 0, 0);
    };

    NeoHooke problem = NeoHooke("", r, C, h, num_cells);
    problem.setup();
    problem.solve();
    return 0;
}
