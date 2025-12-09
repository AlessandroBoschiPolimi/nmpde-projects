#include "NeoHooke.hpp"

#include <deal.II/grid/grid_generator.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/base/tensor.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#define GAMBA_DEBUG false

using namespace pde;

void NeoHooke::setup() {
    std::cout << "===============================" << std::endl;
    
    // Create the mesh
    {
	// Create the mesh using dealii generator, this gives us numbering of the faces
	/* From the documentation:
		
	Faces (quads in 3d): first the two faces with normals in x-, then y- and z-direction. For each two faces: first the face with normal in negative coordinate direction, then the one with normal in positive direction, i.e. the faces are ordered according to their normals pointing in -x, x, -y, y, -z, z direction.
    
	Therefore, the faces are numbered in the ordering: left, right, front, back, bottom and top face:
	
	*       *-------*        *-------*
	*      /|       |       /       /|
	*     / |   3   |      /   5   / |
	*    /  |       |     /       /  |
	*   *   |       |    *-------*   |
	*   | 0 *-------*    |       | 1 *
	*   |  /       /     |       |  /
	*   | /   4   /      |   2   | /
	*   |/       /       |       |/
	*   *-------*        *-------*
	* 
	*/
	std::cout << "Initializing the mesh" << std::endl;
	GridGenerator::subdivided_hyper_cube(mesh, num_cells, 0.0, 1.0, true);
	std::cout << "Number of elements = " << mesh.n_active_cells()
		  << std::endl;
    }

    std::cout << "------------------------------------" << std::endl;

    // Initialize the finite element space
    {
	// TODO: Check that this Finite Element Space suffices
	fe = std::make_unique<FESystem<dim>>(FE_Q<dim>(r)^dim);

	std::cout << "  Degree                     = " << fe->degree << std::endl;
	std::cout << "  DoFs per cell              = " << fe->dofs_per_cell
		  << std::endl;

	// TODO: Check that these quadrature are correct enough
	quadrature = std::make_unique<QGaussSimplex<dim>>(r + 1);
	quadrature_boundary = std::make_unique<QGaussSimplex<dim - 1>>(r + 1);

	std::cout << "  Quadrature points per cell = " << quadrature->size()
	  << std::endl;
    }

    std::cout << "----------------------------------" << std::endl;

    // Initialize the DoF handler.
    {
	std::cout << "Initializing the DoF handler" << std::endl;

	// Initialize the DoF handler with the mesh we constructed.
	dof_handler.reinit(mesh);

	// "Distribute" the degrees of freedom. For a given finite element space,
	// initializes info on the control variables (how many they are, where
	// they are collocated, their "global indices", ...).
	dof_handler.distribute_dofs(*fe);

	std::cout << "  Number of DoFs = " << dof_handler.n_dofs() << std::endl;
    }

    std::cout << "-----------------------------------------------" << std::endl;

    // Initialize the linear system.
    {
	std::cout << "Initializing the linear system" << std::endl;

	std::cout << "  Initializing the sparsity pattern" << std::endl;
	DynamicSparsityPattern dsp(dof_handler.n_dofs());
	DoFTools::make_sparsity_pattern(dof_handler, dsp);
	sparsity_pattern.copy_from(dsp);

	// Then, we use the sparsity pattern to initialize the system matrix
	std::cout << "  Initializing the system matrix" << std::endl;
	jacobian_matrix.reinit(sparsity_pattern);

	// Finally, we initialize the right-hand side and solution vectors.
	std::cout << "  Initializing the system right-hand side" << std::endl;
	residual_vector.reinit(dof_handler.n_dofs());
	std::cout << "  Initializing the solution vector" << std::endl;
	solution.reinit(dof_handler.n_dofs());
	delta.reinit(dof_handler.n_dofs());
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
    residual_vector 	  = 0.0;

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
		std::cout << "fe_values[displacement].get_function_values " << 
		    solution_loc[q] << std::endl;
		std::cout << "fe_values[displacement].get_function_gradients " << 
		   solution_gradient_loc[q] << std::endl;
	}
	#endif

	for ( const unsigned int q : fe_values.quadrature_point_indices() ) {
	    for ( const unsigned int i : fe_values.dof_indices() ) {
		// base to describe v
		const Tensor<2,dim> phi_i_grad = fe_values[displacement].gradient(i, q);

		Tensor<2,dim> displacement_tensor({{1,0,0},{0,1,0},{0,0,1}});
		displacement_tensor += solution_gradient_loc[q];
		const Tensor<2,dim> inverse_displacement = transpose(invert(displacement_tensor));

		for ( const unsigned j : fe_values.dof_indices() ) {
		    // base to describe delta
		    const Tensor<2, dim> phi_j_grad = fe_values[displacement].gradient(j,q);
		    const Tensor<2,dim> second_member = inverse_displacement * transpose(phi_j_grad) * inverse_displacement;

		    // TODO: check correctness of the following multiplication
		    // and if there is an overloaded operator to do it
		    cell_matrix(i,j) += mu * (
			double_contract<0,0,1,1>(phi_j_grad, phi_i_grad) + 
			double_contract<0,0,1,1>(second_member , phi_i_grad)
			) * fe_values.JxW(q);
		    

		}

		cell_rhs(i) -= mu * (
		    double_contract<0,0,1,1>(displacement_tensor, phi_i_grad) -
		    double_contract<0,0,1,1>(inverse_displacement, phi_i_grad)
		) * fe_values.JxW(q);
	    }
	}
	
	if(cell->at_boundary()) {
	    for (unsigned int face_number = 0; 
		face_number < cell->n_faces();
		face_number++) {

		const unsigned int bound_id = cell->face(face_number)->boundary_id();

		if (!cell->face(face_number)->at_boundary()
			|| bound_id == 4
			|| bound_id == 5)
			continue;

		// If current face lies on the boundary, and its boundary ID (or
		// tag) is that of one of the Neumann boundaries, we assemble the
		// boundary integral.

		fe_values_boundary.reinit(cell, face_number);

		for (unsigned int q : fe_values_boundary.quadrature_point_indices()) {
		    for (unsigned int i = 0; i < dofs_per_cell; ++i) {

			#if GAMBA_DEBUG
			    std::cout << "Point(q) " << fe_values_boundary.quadrature_point(q)  << std::endl;
			    std::cout << "fe_values_boundary[displacement].value(i, q) " << fe_values_boundary[displacement].value(i, q) << std::endl;
			#endif

			cell_rhs(i) +=
			h(fe_values_boundary.quadrature_point(q)) *
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
    
    #if GAMBA_DEBUG
	std::cout << "system_rhs " << residual_vector << std::endl;
    #endif

    // TODO: add dirichelet boundary conditions

    //    {
	// std::map<types::global_dof_index, double> boundary_values;
	//
	// std::map<types::boundary_id, const Function<dim> *> boundary_functions;
	// Functions::ZeroFunction<dim>                        zero_function;
	//
	// for (unsigned int i = 0; i < 6; ++i)
	//     boundary_functions[i] = &zero_function;
	//
	// VectorTools::interpolate_boundary_values(dof_handler,
	// 				 boundary_functions,
	// 				 boundary_values);
	//
	// MatrixTools::apply_boundary_values(boundary_values, system_matrix, delta, system_rhs, true);
	//    }

}

void NeoHooke::solve_system() {
    // TODO: Change number of iterations
    SolverControl solver_control(2000, 1e-6 * residual_vector.l2_norm());

    SolverGMRES<Vector<double>> solver(solver_control);
    
    PreconditionSSOR preconditioner;
    preconditioner.initialize(jacobian_matrix);

    solver.solve(jacobian_matrix, delta, residual_vector, preconditioner);
    std::cout << "   " << solver_control.last_step() << " GMRES iterations"
    << std::endl;

}

void NeoHooke::solve() {
    std::cout << "===============================================" << std::endl;

    const unsigned int n_max_iters        = 1000;
    const double       residual_tolerance = 1e-6;

    unsigned int n_iter        = 0;
    double       residual_norm = residual_tolerance + 1;

    while (n_iter < n_max_iters && residual_norm > residual_tolerance)
    {
	assemble_system();
	residual_norm = residual_vector.l2_norm();

	std::cout << "Newton iteration " << n_iter << "/" << n_max_iters
	<< " - ||r|| = " << std::scientific << std::setprecision(6)
	<< residual_norm << std::flush;

	// We actually solve the system only if the residual is larger than the
	// tolerance.
	if (residual_norm <= residual_tolerance) break;

	solve_system();
	solution += delta;

	++n_iter;
    }

    std::cout << " < tolerance" << std::endl;
    std::cout << "===============================================" << std::endl;
}

void NeoHooke::output() const {

}
