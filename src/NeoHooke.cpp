#include "NeoHooke.hpp"

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

// ------------ CPP HEADERS ----------
#include <filesystem>
#include <fstream>

#define GAMBA_DEBUG false

using namespace pde;

void NeoHooke::setup() {

    // TODO: Read from mesh
    pcout << "===============================" << std::endl;
    
    // Create the mesh
	{
	    pcout << "Initializing the mesh" << std::endl;
	
		mesh_generator->Generate(mesh);
		
		pcout << "  Number of elements = " << mesh.n_global_active_cells() << std::endl;
	}

    pcout << "------------------------------------" << std::endl;

    // Initialize the finite element space
    {
		// TODO: Check that this Finite Element Space suffices
		pcout << "Initializing the finite element space" << std::endl;

		if (mesh_generator->ElementType() == Type::Tetrahedra)
			fe = std::make_unique<FESystem<dim>>(FE_SimplexP<dim>(r)^dim);
		else
			fe = std::make_unique<FESystem<dim>>(FE_Q<dim>(r)^dim);

		pcout << "  Degree                     = " << fe->degree << std::endl;
		pcout << "  DoFs per cell              = " << fe->dofs_per_cell
			<< std::endl;

		// TODO: Check that these quadrature are correct enough
		quadrature = std::make_unique<QGauss<dim>>(r + 1);
		quadrature_boundary = std::make_unique<QGauss<dim - 1>>(r + 1);

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
    const unsigned int dofs_per_cell	= fe->dofs_per_cell;
    const unsigned int n_q		= quadrature->size();

    FEValues<dim> fe_values(
		*fe,
		*quadrature,
		update_values | update_gradients |
		update_quadrature_points | update_JxW_values
    );

    FEFaceValues<dim> fe_values_boundary(
		*fe,
		*quadrature_boundary,
		update_values | update_quadrature_points |
		update_JxW_values
    );


    FullMatrix<double>	cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>	cell_rhs(dofs_per_cell);

    std::vector<types::global_cell_index> dof_indices(dofs_per_cell);

    jacobian_matrix = 0.0;
    residual_vector = 0.0;

    std::vector<Tensor<1,dim, double>>	solution_loc(n_q);
    std::vector<Tensor<2,dim, double>>	solution_gradient_loc(n_q);
    

    // Getting a view of the <1,dim> tensor from 0:dim-1
    // An FEValuesExtractor defines a view of the vector
    // (in this case the view contains the exact values we are using
    // due to only needing displacement)
    // on which then the FEValues can perform operations such as
    // in this case to compute gradients of vectors of vectors
    const FEValuesExtractors::Vector displacement(0);


    for ( const auto &cell : dof_handler.active_cell_iterators() ) {
		fe_values.reinit(cell);
		cell_matrix = 0.0;
		cell_rhs    = 0.0;

		fe_values[displacement].get_function_values(solution, solution_loc);
		fe_values[displacement].get_function_gradients(solution, solution_gradient_loc);

		#if GAMBA_DEBUG
		for ( const unsigned int q : fe_values.quadrature_point_indices() ) {
			pcout << "fe_values[displacement].get_function_values " << 
				solution_loc[q] << std::endl;
			pcout << "fe_values[displacement].get_function_gradients " << 
			solution_gradient_loc[q] << std::endl;
		}
		#endif

		for ( const unsigned int q : fe_values.quadrature_point_indices() ) {
			// Compute the displacement tensor I + grad(d(k)) at quadr point q
			Tensor<2,dim> displacement_tensor({{1,0,0},{0,1,0},{0,0,1}});
			displacement_tensor += solution_gradient_loc[q];
			const Tensor<2,dim> inverse_transpose_displacement = transpose(invert(displacement_tensor));
			//Compute determinant of the displacement tensor F
			//I am not sure about the absolute value here, but since we put it into a log, we can get a NaN easily otherwise
			const double determinant_displacement = std::abs(determinant(displacement_tensor));

			for ( const unsigned int i : fe_values.dof_indices() ) {
			// base to describe v
			const Tensor<2,dim> phi_i_grad = fe_values[displacement].gradient(i, q);

				for ( const unsigned j : fe_values.dof_indices() ) {
					// base to describe delta
					const Tensor<2, dim> phi_j_grad = fe_values[displacement].gradient(j,q);
					const Tensor<2,dim> second_member = inverse_transpose_displacement * transpose(phi_j_grad) * inverse_transpose_displacement;

					// TODO: check correctness of the following multiplication
					// and if there is an overloaded operator to do it
					cell_matrix(i,j) += mu * (
							double_contract<0,0,1,1>(phi_j_grad, phi_i_grad) + 
							double_contract<0,0,1,1>(second_member , phi_i_grad)
						) * fe_values.JxW(q);

					//Add two additional terms arrising when lambda =/= 0
					
					cell_matrix(i,j) += lambda * (
							double_contract<0,0,1,1>(phi_j_grad, inverse_transpose_displacement) * 
							double_contract<0,0,1,1>(inverse_transpose_displacement , phi_i_grad)
						) * fe_values.JxW(q);

					//Is the std log relly the best here?
					cell_matrix(i,j) -= lambda * std::log(determinant_displacement) *
						double_contract<0,0,1,1>(second_member,phi_i_grad)
						* fe_values.JxW(q);
				}

			cell_rhs(i) -= mu * (
					double_contract<0,0,1,1>(displacement_tensor, phi_i_grad) -
					double_contract<0,0,1,1>(inverse_transpose_displacement, phi_i_grad)
				) * fe_values.JxW(q);

			//We also have to add the part with lambda to the right side
			cell_rhs(i) -= lambda * std::log(determinant_displacement) * (
					double_contract<0,0,1,1>(inverse_transpose_displacement, phi_i_grad)
				) * fe_values.JxW(q);
			}
		}
		
		if (cell->at_boundary()) {
			for (unsigned int face_number = 0; 
				face_number < cell->n_faces();
				face_number++) {

				const unsigned int bound_id = cell->face(face_number)->boundary_id();

				if (!cell->face(face_number)->at_boundary())
					continue;
				if (!mesh_generator->OnNeumannBoundary(bound_id))
					continue;

				// If current face lies on the boundary, and its boundary ID (or
				// tag) is that of one of the Neumann boundaries, we assemble the
				// boundary integral.

				fe_values_boundary.reinit(cell, face_number);

				for (unsigned int q : fe_values_boundary.quadrature_point_indices()) {
					for (unsigned int i = 0; i < dofs_per_cell; ++i) {
						#if GAMBA_DEBUG
							pcout << "Point(q) " << fe_values_boundary.quadrature_point(q)  << std::endl;
							pcout << "fe_values_boundary[displacement].value(i, q) " << fe_values_boundary[displacement].value(i, q) << std::endl;
						#endif

						cell_rhs(i) +=
							neumann_conds(fe_values_boundary.quadrature_point(q)) *
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
    
    #if GAMBA_DEBUG
	pcout << "system_rhs " << residual_vector << std::endl;
    #endif


    {
		std::map<types::global_dof_index, double> boundary_values;

		VectorTools::interpolate_boundary_values(dof_handler,
					dirichelet_conds,
					boundary_values);

		MatrixTools::apply_boundary_values(boundary_values, jacobian_matrix, delta_owned, residual_vector, false);
    }
}

void NeoHooke::solve_system() {

    SolverControl solver_control(1000, 1e-6 * residual_vector.l2_norm());

    SolverGMRES<TrilinosWrappers::MPI::Vector> solver(solver_control);
    
    TrilinosWrappers::PreconditionSOR preconditioner;
    preconditioner.initialize(jacobian_matrix);

    solver.solve(jacobian_matrix, delta_owned, residual_vector, preconditioner);
    pcout << "   " << solver_control.last_step() << " GMRES iterations"
    	<< std::endl;
}

void NeoHooke::solve() {
    pcout << "===============================================" << std::endl;

    const unsigned int n_max_iters        = 1000;
    const double       residual_tolerance = 1e-6;

    unsigned int n_iter        = 0;
    double       residual_norm = residual_tolerance + 1;

    while (n_iter < n_max_iters && residual_norm > residual_tolerance)
    {
		assemble_system();
		residual_norm = residual_vector.l2_norm();

		pcout << "Newton iteration " << n_iter << "/" << n_max_iters
		<< " - ||r|| = " << std::scientific << std::setprecision(6)
		<< residual_norm << std::flush;

		// We actually solve the system only if the residual is larger than the
		// tolerance.
		if (residual_norm <= residual_tolerance) break;

		solve_system();
		solution_owned += delta_owned;
		solution = solution_owned;

		++n_iter;
    }

    pcout << " < tolerance" << std::endl;
    pcout << "===============================================" << std::endl;
}

void NeoHooke::output() const {
    pcout << "===============================================" << std::endl;

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

    // Use std::filesystem to construct the output file name based on the
    // mesh file name.
    const std::string output_file_name = "../out/output";

    // Finally, we need to write in a format that supports parallel output. This
    // can be achieved in multiple ways (e.g. XDMF/H5). We choose VTU/PVTU files.
    data_out.write_vtu_with_pvtu_record(/* folder = */ "./",
				      /* basename = */ output_file_name,
				      /* index = */ 0,
				      MPI_COMM_WORLD);

    pcout << "Output written to " << output_file_name << std::endl;

    pcout << "===============================================" << std::endl;
}
