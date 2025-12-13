#include "NeoHooke.hpp"
#include "Guccione.hpp"
#include "TestConditions.hpp"


const pde::ForcingTermType select_forcing_term(std::string boolean_value) {
	using namespace pde::TestForcingFunctions;
	return std::stoi(boolean_value) ? bend_rod : null_forcing_term;
}

int main(int argc, char *argv[])
{
    using namespace pde;
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
	
    // Initializing arguments
    const std::string mesh_choice = std::string(argv[1]);
    const std::string neumann_func = std::string(argv[2]);
    const std::string forcing_func = std::string(argv[3]);

    const unsigned int r = 1;
    //Material constant
    const double C = 1; //Pa
    const double lambda = 2; //Pa

    // ---------- INITIALIZE NEUMANN FUNCTION --------------
    TestNeumannConditions::initialize();
    // ----------- CHOOSE NEUMANN FUNCTION ----------------
    const auto h  = TestNeumannConditions::choose_neumann_function(neumann_func);


    std::map<types::boundary_id, const Function<dim> *> boundary_functions;
    Functions::ZeroFunction<dim> zero_function(dim);

    const auto forcing_term = select_forcing_term(forcing_func);

    std::unique_ptr<MechanicalDisplacement> problem;

    if(mesh_choice == "cube") {
	boundary_functions[4] = &zero_function;
	boundary_functions[5] = &zero_function;
	problem = std::make_unique<NeoHooke>(
	    std::make_unique<CubeGenerator<dim>>(), r, boundary_functions, h, forcing_term, C, lambda
	);
    } else {
	// Setting left and right to be still
	boundary_functions[RodGenerator<dim>::left_id] = &zero_function;
	boundary_functions[RodGenerator<dim>::right_id] = &zero_function;

	problem = std::make_unique<NeoHooke>(
	    std::make_unique<RodGenerator<dim>>(), r, boundary_functions, h, forcing_term, C, lambda
	);
    }

    problem->setup();
    problem->solve();
    problem->output();
    return 0;
}
