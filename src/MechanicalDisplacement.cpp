#include "headers/MechanicalDisplacement.hpp"

namespace pde {

void MechanicalDisplacement::setup() {
    pcout << "===============================" << std::endl;
    
    // Create the mesh
	{
	    pcout << "Initializing the mesh" << std::endl;
	
		config.mesh_generator->Generate(mesh);
		
		pcout << "  Number of elements = " << mesh.n_global_active_cells() << std::endl;
	}

    pcout << "------------------------------------" << std::endl;

    // Initialize the finite element space
    {
		// TODO: Check that this Finite Element Space suffices
		pcout << "Initializing the finite element space" << std::endl;

		if (config.mesh_generator->ElementType() == Type::Tetrahedra)
			fe = std::make_unique<FESystem<dim>>(FE_SimplexP<dim>(config.r)^dim);
		else
			fe = std::make_unique<FESystem<dim>>(FE_Q<dim>(config.r)^dim);

		pcout << "  Degree                     = " << fe->degree << std::endl;
		pcout << "  DoFs per cell              = " << fe->dofs_per_cell << std::endl;

		// TODO: Check that these quadrature are correct enough
		if (config.mesh_generator->ElementType() == Type::Tetrahedra) {
			quadrature = std::make_unique<QGaussSimplex<dim>>(config.r + 1);
			quadrature_boundary = std::make_unique<QGaussSimplex<dim - 1>>(config.r + 1);
		} else {
			quadrature = std::make_unique<QGauss<dim>>(config.r + 1);
			quadrature_boundary = std::make_unique<QGauss<dim - 1>>(config.r + 1);
		}

		pcout << "  Quadrature points per cell = " << quadrature->size() << std::endl;
    }

    pcout << "----------------------------------" << std::endl;

    // Initialize the DoF handler.
    {
		pcout << "Initializing the DoF handler" << std::endl;

		// Initialize the DoF handler with the mesh we constructed.
		dof_handler.reinit(mesh);

		// "Distribute" the degrees of freedom. For a given finite element space,
		// initializes info on the control variables (how many they are, where
		// they are collocated, their "global indices", ...).
		dof_handler.distribute_dofs(*fe);

		// We retrieve the set of locally owned DoFs, which will be useful when
		// initializing linear algebra classes.
		locally_owned_dofs = dof_handler.locally_owned_dofs();
		locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

		pcout << "  Number of DoFs = " << dof_handler.n_dofs() << std::endl;
    }

    pcout << "-----------------------------------------------" << std::endl;

    // Initialize the linear system.
    {
		pcout << "Initializing the linear system" << std::endl;

		pcout << "  Initializing the sparsity pattern" << std::endl;

		// To initialize the sparsity pattern, we use Trilinos' class, that manages
		// some of the inter-process communication.
		TrilinosWrappers::SparsityPattern sparsity_pattern(locally_owned_dofs, MPI_COMM_WORLD);
		DoFTools::make_sparsity_pattern(dof_handler, sparsity_pattern);

		// After initialization, we need to call compress, so that all process
		// retrieve the information they need for the rows they own (i.e. the rows
		// corresponding to locally owned DoFs).
		sparsity_pattern.compress();

		// Then, we use the sparsity pattern to initialize the system matrix
		pcout << "  Initializing the system matrix" << std::endl;
		jacobian_matrix.reinit(sparsity_pattern);

		// Finally, we initialize the right-hand side and solution vectors.
		pcout << "  Initializing the system right-hand side" << std::endl;
		residual_vector.reinit(locally_owned_dofs, MPI_COMM_WORLD);
		pcout << "  Initializing the solution vector" << std::endl;
		solution_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);
		solution.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
		delta_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    }
}


void MechanicalDisplacement::solve() {
    pcout << "===============================================" << std::endl;

    const unsigned int n_max_iters        = 1000;
    const double       residual_tolerance = 1e-6;

    unsigned int n_iter        = 0;
    double       residual_norm = residual_tolerance + 1;

	// apply the required D conditions on the initial guess
	// apply_dirchlet_to_initial_solution();

    while (n_iter < n_max_iters && residual_norm > residual_tolerance)
    {
		assemble_system();
		residual_norm = residual_vector.l2_norm();

		pcout << "Newton iteration " << n_iter << "/" << n_max_iters
			  << " - ||r|| = " << std::scientific << std::setprecision(6)
			  << residual_norm << std::flush;

		// We actually solve the system only if the residual is larger than the tolerance.
		if (residual_norm <= residual_tolerance)
			break;
		
		// apply a homogeneous D condition on the guess update, to maintain the correct boundary condition.
		// apply_zero_dirchlet_to_newton_update();

		solve_system();
		
		if (config.newton_damping)
			delta_owned *= config.newton_scaling;
		solution_owned += delta_owned;
		solution = solution_owned;
		solution.update_ghost_values();

		++n_iter;
		// output(n_iter);
    }

    pcout << " < tolerance" << std::endl;
    pcout << "===============================================" << std::endl;
}


void MechanicalDisplacement::apply_dirchlet_to_initial_solution()
{
	std::map<types::global_dof_index, double> boundary_values;
	VectorTools::interpolate_boundary_values(dof_handler, config.dirichelet_conds, boundary_values);
	for (const auto &bc : boundary_values)
		if (solution_owned.locally_owned_elements().is_element(bc.first))
			solution_owned[bc.first] = bc.second;
	solution = solution_owned;
	solution.update_ghost_values();

	// VectorTools::interpolate_boundary_values(dof_handler, dirichelet_conds, solution);

	// assemble_system();

	// VectorTools::interpolate_boundary_values(dof_handler, dirichelet_conds, boundary_values);
	// MatrixTools::apply_boundary_values(boundary_values, jacobian_matrix, solution, residual_vector);
}

void MechanicalDisplacement::apply_zero_dirchlet_to_newton_update()
{
	std::map<types::global_dof_index, double> boundary_values;
	std::map<types::boundary_id, const Function<dim> *> zero_dirichelet_conds;
	Functions::ZeroFunction<dim> zero(dim);

	for (auto& d : config.dirichelet_conds)
		zero_dirichelet_conds[d.first] = &zero;

	delta_owned = 0.0;
	VectorTools::interpolate_boundary_values(dof_handler, zero_dirichelet_conds, boundary_values);
	MatrixTools::apply_boundary_values(boundary_values, jacobian_matrix, delta_owned, residual_vector, false);
}


void MechanicalDisplacement::solve_system() {
    SolverControl solver_control(config.iterations, 1e-6 * residual_vector.l2_norm());

    SolverGMRES<TrilinosWrappers::MPI::Vector> solver(solver_control);
    
    TrilinosWrappers::PreconditionAMG preconditioner;
    preconditioner.initialize(jacobian_matrix);

    solver.solve(jacobian_matrix, delta_owned, residual_vector, preconditioner);
    pcout << "   " << solver_control.last_step() << " GMRES iterations" << std::endl;
}


void MechanicalDisplacement::output(int ts) const {
    pcout << "===============================================" << std::endl;

	{
		std::filesystem::path p(config.output_filename);
		std::filesystem::create_directories(p.parent_path());
	}
	
    std::vector<std::string> solution_names(dim, "displacement");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation
    (dim, DataComponentInterpretation::component_is_part_of_vector);

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(dof_handler, solution, solution_names, data_component_interpretation);

    // Partitioning data in order to write in parallel
    std::vector<unsigned int> partition_int(mesh.n_active_cells());
    GridTools::get_subdomain_association(mesh, partition_int);
    const Vector<double> partitioning(partition_int.begin(), partition_int.end());
    data_out.add_data_vector(partitioning, "partitioning");

    data_out.build_patches();

    std::filesystem::path output_path = config.output_filename;
	std::filesystem::path parent_path = output_path.parent_path();
	parent_path /= ""; // ensures '/' at the end of parent_path
	std::filesystem::path filename = output_path.filename();

    data_out.write_vtu_with_pvtu_record(parent_path.string(), filename.string(), ts, MPI_COMM_WORLD);

    pcout << "Output written to " << config.output_filename << std::endl;

    pcout << "===============================================" << std::endl;
}

}