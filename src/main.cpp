#include "NeoHooke.hpp"

int main(int argc, char *argv[])
{
    using namespace pde;
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

    const unsigned int dim = MechanicalDisplacement::dim;
    const unsigned int r = 1;

    //Material constant
    const double C = 1; //Pa
    const double lambda = 2; //Pa
    //Number of cells inside a hypercube
    const unsigned int num_cells = 10;
    const double tau_0 = 0.5;
    //Pulling the cube appart <-|  |-> in x,y directions
    //The dirichlet conditions are set to 0 on the z boundaries
    //The object can change volume, since lambda = 0 in our model
    const auto h  = [&tau_0](const Point<dim> &p) {
	double small_tol = 1e-13;
	if (std::abs(p[0]) < small_tol && std::abs(p[1]) > small_tol)
		return Point<dim>(tau_0, 0, 0);
	else if (std::abs(p[0]) > (1 - small_tol) && std::abs(p[1]) > small_tol)
		return Point<dim>(-tau_0, 0, 0);
	else if (std::abs(p[1]) < small_tol && std::abs(p[0]) > small_tol)
		return Point<dim>(0, 0, 0);
	else if (std::abs(p[1]) > (1 - small_tol) && std::abs(p[0]) > small_tol)
		return Point<dim>(0, 0, 0);
	else
		return Point<dim>(0, 0, 0);
    };

    
    std::map<types::boundary_id, const Function<dim> *> boundary_functions;
    Functions::ZeroFunction<dim> zero_function(dim);

    boundary_functions[4] = &zero_function;
    boundary_functions[5] = &zero_function;

    NeoHooke problem = NeoHooke("gds", r, boundary_functions, h, num_cells, C, lambda);
    problem.setup();
    problem.solve();
    problem.output();
    return 0;
}
