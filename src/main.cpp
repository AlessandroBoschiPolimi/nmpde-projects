#include "NeoHooke.hpp"
#include "Guccione.hpp"
#include "TestFunctions.hpp"

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

    // ---------- INITIALIZE NEUMANN FUNCTION --------------
    TestFunctions::initialize();
    // ----------- CHOOSE NEUMANN FUNCTION ----------------
    const auto h  = TestFunctions::choose_neumann_function(argv[2]);


    std::map<types::boundary_id, const Function<dim> *> boundary_functions;
    Functions::ZeroFunction<dim> zero_function(dim);

    
    // TODO: fix this

    std::unique_ptr<MechanicalDisplacement> problem;

    if(std::string(argv[1]) == "cube") {
	boundary_functions[4] = &zero_function;
	boundary_functions[5] = &zero_function;
	problem = std::make_unique<NeoHooke>(
	    std::make_unique<CubeGenerator<dim>>(), r, boundary_functions, h, num_cells, C, lambda
	);
    } else {
	// Setting left and right to be still
	boundary_functions[RodGenerator<dim>::left_id] = &zero_function;
	boundary_functions[RodGenerator<dim>::right_id] = &zero_function;

	problem = std::make_unique<NeoHooke>(
	    std::make_unique<RodGenerator<dim>>(), r, boundary_functions, h, num_cells, C, lambda
	);
    }

    problem->setup();
    problem->solve();
    problem->output();
    return 0;
}
