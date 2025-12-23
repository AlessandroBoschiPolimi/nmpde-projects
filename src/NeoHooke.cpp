#include "headers/NeoHooke.hpp"

// --------- DEALII HEADERS ----------
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>

#include <deal.II/base/tensor.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/numerics/data_out.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/fe/mapping_fe.h>

#include "headers/TestConditions.hpp"


#define GAMBA_DEBUG false

using namespace pde;

void NeoHooke::setup() {
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


void NeoHooke::assemble_system() {
    const unsigned int dofs_per_cell = fe->dofs_per_cell;
    const unsigned int n_q           = quadrature->size();

	// MappingFE<dim> mapping(FE_SimplexP<dim>(2));

    FEValues<dim> fe_values(
		// mapping,
		*fe,
		*quadrature,
		update_values | update_gradients |
		update_quadrature_points | update_JxW_values
    );

    FEFaceValues<dim> fe_values_boundary(
		// mapping,
		*fe,
		*quadrature_boundary,
		update_values | update_quadrature_points | update_JxW_values
    );


    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);

    std::vector<types::global_cell_index> dof_indices(dofs_per_cell);

    jacobian_matrix = 0.0;
    residual_vector = 0.0;

    std::vector<Tensor<1, dim, double>> solution_loc(n_q);
    std::vector<Tensor<2, dim, double>> solution_gradient_loc(n_q);
    

    // Getting a view of the <1,dim> tensor from 0:dim-1
    // An FEValuesExtractor defines a view of the vector
    // (in this case the view contains the exact values we are using
    // due to only needing displacement)
    // on which then the FEValues can perform operations such as
    // in this case to compute gradients of vectors of vectors
    const FEValuesExtractors::Vector displacement(0);


    for (const auto &cell : dof_handler.active_cell_iterators()) {
		if (!cell->is_locally_owned())
			continue;

		fe_values.reinit(cell);
		cell_matrix = 0.0;
		cell_rhs    = 0.0;

		fe_values[displacement].get_function_values(solution, solution_loc);
		fe_values[displacement].get_function_gradients(solution, solution_gradient_loc);

		for (const unsigned int q : fe_values.quadrature_point_indices()) {
			// Compute the displacement tensor I + grad(d(k)) at quadr point q
			Tensor<2, dim> displacement_tensor({{1,0,0},{0,1,0},{0,0,1}});
			displacement_tensor += solution_gradient_loc[q];
			const Tensor<2, dim> inverse_displacement = invert(displacement_tensor);
			const Tensor<2, dim> inverse_transpose_displacement = transpose(inverse_displacement);
			// Compute determinant of the displacement tensor F
			// I am not sure about the absolute value here, but since we put it into a log, we can get a NaN easily otherwise
			const double determinant_displacement = determinant(displacement_tensor);
			Assert(determinant_displacement > 0, ExcMessage("det(F) <= 0"));

			for (const unsigned int i : fe_values.dof_indices()) {
				// base to describe v
				const Tensor<2, dim> phi_i_grad = fe_values[displacement].gradient(i, q);

				for (const unsigned j : fe_values.dof_indices()) {
					// base to describe delta
					const Tensor<2, dim> phi_j_grad = fe_values[displacement].gradient(j, q);

					#if true // chatgpt version

					// trace(F^{-1} dF)
					double tr_Finv_dF = scalar_product(inverse_displacement, phi_j_grad);

					// term 1: mu * dF
					double term1 = mu * scalar_product(phi_j_grad, phi_i_grad);

					// term 2: mu * F^{-T} dF^T F^{-T}
					Tensor<2, dim> A = inverse_transpose_displacement * transpose(phi_j_grad) * inverse_transpose_displacement;
					double term2 = mu * scalar_product(A, phi_i_grad);

					// term 3: lambda * tr(F^{-1} dF) * F^{-T}
					double term3 = lambda * tr_Finv_dF * scalar_product(inverse_transpose_displacement, phi_i_grad);

					// term 4: - lambda * log(J) * F^{-T} dF^T F^{-T}
					double term4 = -lambda * std::log(determinant_displacement) * scalar_product(A, phi_i_grad);

					cell_matrix(i, j) += (term1 + term2 + term3 + term4) * fe_values.JxW(q);

					#else // our version

					const Tensor<2, dim> second_member = inverse_transpose_displacement * transpose(phi_j_grad) * inverse_transpose_displacement;

					// TODO: check correctness of the following multiplication
					// and if there is an overloaded operator to do it
					cell_matrix(i,j) += mu * (
							double_contract<0,0,1,1>(phi_j_grad, phi_i_grad) + 
							double_contract<0,0,1,1>(second_member, phi_i_grad)
						) * fe_values.JxW(q);

					// Add two additional terms arrising when lambda =/= 0
					
					cell_matrix(i,j) += lambda * (
							double_contract<0,0,1,1>(phi_j_grad, inverse_transpose_displacement) * 
							double_contract<0,0,1,1>(inverse_transpose_displacement, phi_i_grad)
						) * fe_values.JxW(q);

					// Is the std log relly the best here?
					cell_matrix(i,j) -= lambda * std::log(determinant_displacement) *
						double_contract<0,0,1,1>(second_member, phi_i_grad)
						* fe_values.JxW(q);

					#endif
				}

				cell_rhs(i) -= mu * (
						double_contract<0,0,1,1>(displacement_tensor, phi_i_grad) -
						double_contract<0,0,1,1>(inverse_transpose_displacement, phi_i_grad)
					) * fe_values.JxW(q);

				// We also have to add the part with lambda to the right side
				cell_rhs(i) -= lambda * std::log(determinant_displacement) * (
						double_contract<0,0,1,1>(inverse_transpose_displacement, phi_i_grad)
					) * fe_values.JxW(q);

				cell_rhs(i) += config.forcing_term(fe_values.quadrature_point(q)) *
						fe_values[displacement].value(i, q)[1] * (int)(fe_values.quadrature_point(q)[0] > 0) * fe_values.JxW(q);
			}
		}
		
		if (cell->at_boundary()) {
			for (unsigned int face_number = 0; face_number < cell->n_faces(); face_number++) {
				const unsigned int bound_id = cell->face(face_number)->boundary_id();

				if (!cell->face(face_number)->at_boundary())
					continue;
				if (!config.neumann_ids.contains(bound_id))
					continue;

				fe_values_boundary.reinit(cell, face_number);

				for (unsigned int q : fe_values_boundary.quadrature_point_indices()) {
					for (unsigned int i = 0; i < dofs_per_cell; ++i) {
						#if GAMBA_DEBUG
							pcout << "Point(q) " << fe_values_boundary.quadrature_point(q)  << std::endl;
							pcout << "Neumann Conditions: " <<
								config.neumann_conds(fe_values_boundary.quadrature_point(q)) << std::endl;
							pcout << "fe_values_boundary[displacement].value(i, q) " << 
								fe_values_boundary[displacement].value(i, q) << std::endl;
						#endif

						cell_rhs(i) +=
							config.neumann_conds(fe_values_boundary.quadrature_point(q)) *
							fe_values_boundary[displacement].value(i, q) *     
							fe_values_boundary.JxW(q);
					}
				}
			}
		}
	
		// At this point the local matrix and vector are constructed: we need
		// to sum them into the global matrix and vector. To this end, we need
		// to retrieve the global indices of the DoFs of current cell.
		cell->get_dof_indices(dof_indices);

		jacobian_matrix.add(dof_indices, cell_matrix);
		residual_vector.add(dof_indices, cell_rhs);
    }

    // Synchronize matrix and residual vector values
    jacobian_matrix.compress(VectorOperation::add);
    residual_vector.compress(VectorOperation::add);

	#if true
		std::map<types::global_dof_index, double> boundary_values;
		VectorTools::interpolate_boundary_values(dof_handler, config.dirichelet_conds, boundary_values);

		// for (auto &bc : boundary_values)
		// 	bc.second -= solution[bc.first];

		delta_owned = 0.0;
		MatrixTools::apply_boundary_values(boundary_values, jacobian_matrix, delta_owned, residual_vector, false);
	#endif
}

void NeoHooke::solve_system() {
    SolverControl solver_control(config.iterations, 1e-6 * residual_vector.l2_norm());

    SolverGMRES<TrilinosWrappers::MPI::Vector> solver(solver_control);
    
    TrilinosWrappers::PreconditionAMG preconditioner;
    preconditioner.initialize(jacobian_matrix);

    solver.solve(jacobian_matrix, delta_owned, residual_vector, preconditioner);
    pcout << "   " << solver_control.last_step() << " GMRES iterations" << std::endl;
}

void NeoHooke::apply_dirchlet_to_initial_solution()
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

void NeoHooke::apply_zero_dirchlet_to_newton_update()
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

void NeoHooke::solve() {
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
    }

    pcout << " < tolerance" << std::endl;
    pcout << "===============================================" << std::endl;
}

void NeoHooke::output() const {
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

    data_out.write_vtu_with_pvtu_record(parent_path.string(), filename.string(), 0, MPI_COMM_WORLD);

    pcout << "Output written to " << config.output_filename << std::endl;

    pcout << "===============================================" << std::endl;
}
