#include <iostream>
#include <filesystem>
#include "MechanicalDisplacement.hpp"
#include <cstdlib>
int main(int argc, char *argv[])
{
	std::cout << "Hello" << '\n';
	Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
	const unsigned int r = 2;
	//Number of cells inside a hypercube
	const unsigned int num_cells = 10;
	const double tau_0 = 0.1;
	const auto h  = [&tau_0](const Point<3> &p) {
	if (std::abs(p[0]) < 1e-8 && std::abs(p[1]) > 1e-8)
		return Point<3>(-tau_0, 0, 0);
	else if (std::abs(p[0]) > 1 - 1e-8 && std::abs(p[1]) > 1e-8)
		return Point<3>(tau_0, 0, 0);
	else if (std::abs(p[1]) < 1e-8 && std::abs(p[0]) > 1e-8)
		return Point<3>(0, -tau_0, 0);
	else if (std::abs(p[1]) > 1 - 1e-8 && std::abs(p[0]) > 1e-8)
		return Point<3>(0, tau_0, 0);
	else
		return Point<3>(0, 0, 0);
	};
	MechanicalDisplacement Mech("Nenene_To_Se_Nedela", r, h, num_cells);
	Mech.setup();
	//Mech.solve_newton();
	//Mech.output();
	return 0;
}
