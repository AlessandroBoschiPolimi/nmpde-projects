#include <iostream>
#include <filesystem>
#include "MechanicalDisplacement.hpp"
#include <cstdlib>
int main(int argc, char *argv[])
{
	std::cout << "Hello" << '\n';
	Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
	const unsigned int r = 1;
	//Material constant
	const double C = 3; //Pa
	//Number of cells inside a hypercube
	const unsigned int num_cells = 8;
	const double tau_0 = 0.5;
	//Pulling the cube appart <-|  |-> in x,y directions
	//The dirichlet conditions are set to 0 on the z boundaries
	//The object can change volume, since lambda = 0 in our model
	const auto h  = [&tau_0](const Point<3> &p) {
		double small_tol = 1e-13;
	if (std::abs(p[0]) < small_tol && std::abs(p[1]) > small_tol)
		return Point<3>(-tau_0, 0, 0);
	else if (std::abs(p[0]) > (1 - small_tol) && std::abs(p[1]) > small_tol)
		return Point<3>(tau_0, 0, 0);
	else if (std::abs(p[1]) < small_tol && std::abs(p[0]) > small_tol)
		return Point<3>(0, -tau_0, 0);
	else if (std::abs(p[1]) > (1 - small_tol) && std::abs(p[0]) > small_tol)
		return Point<3>(0, tau_0, 0);
	else
		return Point<3>(0, 0, 0);
	};
	MechanicalDisplacement Mech("Nenene_To_Se_Nedela", r, C, h, num_cells);
	Mech.setup();
	Mech.solve_newton();
	Mech.output();
	return 0;
}
