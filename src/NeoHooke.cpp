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
	Triangulation<dim> mesh_serial;
	std::cout << "Initializing the mesh" << std::endl;
	GridGenerator::subdivided_hyper_cube(mesh, num_cells, 0.0, 1.0, true);
	std::cout << "  Number of elements = " << mesh.n_active_cells()
		  << std::endl;
    }

    std::cout << "------------------------------------" << std::endl;

    // Initialize the finite element space
    {
	fe = std::make_unique<FESystem<dim>>(FE_Q<dim>(r)^dim);

	std::cout << "  Degree                     = " << fe->degree << std::endl;
	std::cout << "  DoFs per cell              = " << fe->dofs_per_cell
		  << std::endl;

	quadrature = std::make_unique<QGaussSimplex<dim>>(r + 1);

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
	system_matrix.reinit(sparsity_pattern);

	// Finally, we initialize the right-hand side and solution vectors.
	std::cout << "  Initializing the system right-hand side" << std::endl;
	system_rhs.reinit(dof_handler.n_dofs());
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

    FullMatrix<double>	cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>	cell_rhs(dofs_per_cell);

    std::vector<types::global_cell_index> dof_indices(dofs_per_cell);

    system_matrix = 0.0;
    system_rhs 	  = 0.0;

    std::vector<Tensor<1,dim, double>>	solution_loc(n_q);
    std::vector<Tensor<2,dim, double>>	solution_gradient_loc(n_q);

    // Getting a view of the <2,dim> tensor from 0:dim-1
    const FEValuesExtractors::Vector displacement(0);
    double value;
    for ( const auto &cell : dof_handler.active_cell_iterators() ) {

	fe_values.reinit(cell);
	cell_matrix = 0.0;
	cell_rhs    = 0.0;

	fe_values[displacement].get_function_values(solution, solution_loc);
	fe_values[displacement].get_function_gradients(solution, solution_gradient_loc);

	for ( const unsigned int i : fe_values.dof_indices() ) {
	    const unsigned int comp_i =
		fe->system_to_component_index(i).first;
	    for ( const unsigned int q : fe_values.quadrature_point_indices() ) {
		Tensor<2,dim> displacement_tensor({{1,0,0},{0,1,0},{0,0,1}});
		displacement_tensor += solution_gradient_loc[q];
		Tensor<2,dim> inverse_displacement = transpose(invert(displacement_tensor));
		for ( const unsigned j : fe_values.dof_indices() ) {
		    const unsigned int comp_j =
			fe->system_to_component_index(j).first;
		    value = 0;
		    for(unsigned int a = 0; a < dim; a++)
			value += inverse_displacement[a][comp_j] * fe_values.shape_grad(i, q)[comp_j] * inverse_displacement[comp_i][a];
		   

		    cell_matrix(i,j) += comp_i == comp_j
			? fe_values.shape_grad(i,q) * fe_values.shape_grad(j,q)
			* fe_values.JxW(q) : 0;

		    cell_matrix(i,j) += value * fe_values.shape_grad(j,q)[comp_i] * fe_values.JxW(q);
		}

		for(unsigned int j = 0; j < dim; j++) {
		    // TODO: Check this computation cause I'm not conviced it is right in any way
		    cell_rhs(i) += displacement_tensor[j][comp_i] * fe_values.shape_grad(i,q)[comp_i] * fe_values.JxW(q);
		    cell_rhs(i) += inverse_displacement[j][comp_i] * fe_values.shape_grad(i,q)[comp_i] * fe_values.JxW(q);
		}

	    }
	}
	
	// At this point the local matrix and vector are constructed: we need
	// to sum them into the global matrix and vector. To this end, we need
	// to retrieve the global indices of the DoFs of current cell.
	cell->get_dof_indices(dof_indices);

	system_matrix.add(dof_indices, cell_matrix);
	system_rhs.add(dof_indices, cell_rhs);
    }

    // TODO: add boundary conditions
    {
	std::map<types::global_dof_index, double> boundary_values;

	std::map<types::boundary_id, const Function<dim> *> boundary_functions;
	Functions::ZeroFunction<dim>                        zero_function;

	for (unsigned int i = 0; i < 6; ++i)
	    boundary_functions[i] = &zero_function;

	VectorTools::interpolate_boundary_values(dof_handler,
					 boundary_functions,
					 boundary_values);

	MatrixTools::apply_boundary_values(boundary_values, system_matrix, delta, system_rhs, true);
    }

}

void NeoHooke::solve_system() {
    SolverControl solver_control(1000, 1e-6 * system_rhs.l2_norm());

    SolverGMRES<Vector<double>> solver(solver_control);
    
    PreconditionSSOR preconditioner;
    preconditioner.initialize(system_matrix);

    solver.solve(system_matrix, delta, system_rhs, preconditioner);
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
	residual_norm = system_rhs.l2_norm();

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
