#include "MechanicalDisplacement.hpp"

#include <deal.II/fe/fe_q.h>

void
MechanicalDisplacement::setup()
{
	// Create the mesh.
	{
		pcout << "Initializing the mesh" << std::endl;

		// First we read the mesh from file into a serial (i.e. not parallel)
		// triangulation.
		Triangulation<dim> mesh_serial;

		{
			GridIn<dim> grid_in;
			grid_in.attach_triangulation(mesh_serial);

			std::ifstream grid_in_file(mesh_file_name);
			grid_in.read_msh(grid_in_file);
		}

		// Then, we copy the triangulation into the parallel one.
		{
			GridTools::partition_triangulation(mpi_size, mesh_serial);
			const auto construction_data = TriangulationDescription::Utilities::create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
			mesh.create_triangulation(construction_data);
		}

		// Notice that we write here the number of *global* active cells (across all
		// processes).
		pcout << "  Number of elements = " << mesh.n_global_active_cells() << std::endl;
	}

	pcout << "-----------------------------------------------" << std::endl;

	// Initialize the finite element space. This is the same as in serial codes.
	{
		pcout << "Initializing the finite element space" << std::endl;

		fe = std::make_unique<FESystem<dim>>(FE_Q<dim>(r)^dim);

		pcout << "  Degree                     = " << fe->degree << std::endl;
		pcout << "  DoFs per cell              = " << fe->dofs_per_cell << std::endl;

		quadrature = std::make_unique<QGaussSimplex<dim>>(r + 1);

		pcout << "  Quadrature points per cell = " << quadrature->size() << std::endl;
	}

	pcout << "-----------------------------------------------" << std::endl;

	// Initialize the DoF handler.
	{
		pcout << "Initializing the DoF handler" << std::endl;

		dof_handler.reinit(mesh);
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
		TrilinosWrappers::SparsityPattern sparsity(locally_owned_dofs, MPI_COMM_WORLD);
		DoFTools::make_sparsity_pattern(dof_handler, sparsity);

		// After initialization, we need to call compress, so that all process
		// retrieve the information they need for the rows they own (i.e. the rows
		// corresponding to locally owned DoFs).
		sparsity.compress();

		// Then, we use the sparsity pattern to initialize the system matrix. Since
		// the sparsity pattern is partitioned by row, so will the matrix.
		pcout << "  Initializing the system matrix" << std::endl;
		jacobian_matrix.reinit(sparsity);

		// Finally, we initialize the right-hand side and solution vectors.
		pcout << "  Initializing the system right-hand side" << std::endl;
		residual_vector.reinit(locally_owned_dofs, MPI_COMM_WORLD);
		pcout << "  Initializing the solution vector" << std::endl;
		solution_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);
		solution.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
		delta_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);
	}
}


void
MechanicalDisplacement::assemble_system()
{
	const unsigned int dofs_per_cell = fe->dofs_per_cell;
	const unsigned int n_q           = quadrature->size();

	FEValues<dim> fe_values(*fe, *quadrature, update_values | update_gradients | update_quadrature_points | update_JxW_values);

	// Since we need to compute integrals on the boundary for Neumann conditions,
	// we also need a FEValues object to compute quantities on boundary edges (faces).
	FEFaceValues<dim> fe_values_boundary(*fe, *quadrature_boundary, update_values | update_quadrature_points | update_JxW_values);

	FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
	/// \int_{\Omega} \sum_j^{dim} P_{\textrm{comp}(i),j} (d^{(k)}) \partial_j [\psi_i]_{\textrm{comp}(i)}
	/// - \int_{\Gamma_N} h_{\textrm{comp}(i)} \psi_i
	Vector<double>     cell_rhs(dofs_per_cell);

	std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

	jacobian_matrix = 0.0;
	residual_vector = 0.0;

	// We use these vectors to store the old solution (i.e. at previous Newton
	// iteration) and its gradient on quadrature nodes of the current cell.
	std::vector<double>         solution_loc(n_q);
	std::vector<Tensor<2, dim>> solution_gradient_loc(n_q);
	FEValuesExtractors::Vector u(0);

	for (const auto &cell : dof_handler.active_cell_iterators())
	{
		if (!cell->is_locally_owned())
			continue;

		fe_values.reinit(cell);

		cell_matrix = 0.0;
		cell_rhs    = 0.0;

		cell->get_dof_indices(dof_indices);

		fe_values.get_function_values(solution, solution_loc);
		// ChatGPT gave me this, i have no idea what is going on, what is a FEValuesExtractor?
		// is there a cleaner way to achieve this? why does the "standard" get_function_gradients() without extractor work with Tensor<1, dim>?
		// why do i have to specify an extractor without any initialization, can't it be done automatically? what am i doing awake at 3:50am?
		Vector<double> solution_at_dofs(dofs_per_cell);
		for (unsigned int i = 0; i < dofs_per_cell; ++i)
			solution_at_dofs[i] = solution[dof_indices[i]];
		fe_values[u].get_function_gradients(solution_at_dofs, solution_gradient_loc);

		const double mu = 1.0;

		for (unsigned int q = 0; q < n_q; ++q)
		{
			// Building P_{nh}
			Tensor<2, dim> I = Tensor<2, dim>();
			for (unsigned int i = 0; i < dim; ++i)
				I[i][i] = 1.0;
			
			auto idk = I - solution_gradient_loc[q];
			auto idk_tinv = invert(transpose(idk));
			
			Tensor<2,dim> piola_m2_politecnico_universita_degli_studi = mu * (idk - idk_tinv);

			for(const unsigned int i : fe_values.dof_indices())
			{
				const unsigned int component_i = fe->system_to_component_index(i).first;
				
				for (const unsigned int j : fe_values.dof_indices())
				{
					const unsigned int component_j = fe->system_to_component_index(j).first;

					// TODO
					cell_matrix(i,j) += (component_i == component_j)
										? fe_values.shape_grad(i, q) * fe_values.shape_grad(j,q) * fe_values.JxW(q)
										: 0;
					
					cell_matrix(i,j) -= (component_i == component_j)
										? mu * fe_values.shape_grad(i, q) * fe_values.shape_grad(j, q) * fe_values.JxW(q)
										: 0;
				}
					

				for (unsigned int j = 0; j < dim; j++)
				{
					/// \int_{\Omega} \sum_j^{dim} P_{\textrm{comp}(i),j} (d^{(k)}) \partial_j [\psi_i]_{\textrm{comp}(i)}
					cell_rhs(i) += piola_m2_politecnico_universita_degli_studi[component_i][j]
								 * fe_values.shape_grad(i, q)[j] // assuming really hard that this returns a Tensor<1, dim> representing the 3 derivatives of the only non zero entry in the base.
								 * fe_values.JxW(q);
				}
			}
		}

		// Temporary function R3 -> R3 identically equal to 1
		Functions::ConstantFunction<dim> h(1., 3);

		// Consider the contrubution of the Neumann boundary conditions on the RHS
		if (cell->at_boundary())
		{
			for (unsigned int face_number = 0; face_number < cell->n_faces(); ++face_number)
			{
				// TODO: what about those ids?
				if (!(cell->face(face_number)->at_boundary() && cell->face(face_number)->boundary_id() == 2))
					continue;

				fe_values_boundary.reinit(cell, face_number);

				for (unsigned int q = 0; q < quadrature_boundary->size(); ++q)
				{
					for (unsigned int i = 0; i < dofs_per_cell; ++i)
					{
						const unsigned int component_i = fe->system_to_component_index(i).first;
					
						Vector<double> hv(dim);
						h.vector_value(fe_values_boundary.quadrature_point(q), hv);
						
						/// -\int_{\Gamma_N} h_{\textrm{comp}(i)} \psi_i
						cell_rhs(i) -= hv[component_i]
									 * fe_values_boundary.shape_value(i, q)
									 * fe_values_boundary.JxW(q);
						
					}
				}
			}
		}

		cell->get_dof_indices(dof_indices);

		jacobian_matrix.add(dof_indices, cell_matrix);
		residual_vector.add(dof_indices, cell_rhs);
	}

	jacobian_matrix.compress(VectorOperation::add);
	residual_vector.compress(VectorOperation::add);

	// TODO: Boundary conditions, now it's not zero on the Ditto boundary, but an arbitrary g.
	{
		std::map<types::global_dof_index, double> boundary_values;

		std::map<types::boundary_id, const Function<dim> *> boundary_functions;
		Functions::ZeroFunction<dim>                        zero_function;

		for (unsigned int i = 0; i < 6; ++i)
			boundary_functions[i] = &zero_function;

		VectorTools::interpolate_boundary_values(dof_handler, boundary_functions, boundary_values);
		MatrixTools::apply_boundary_values(boundary_values, jacobian_matrix, delta_owned, residual_vector, true);
	}
}

void
MechanicalDisplacement::solve_system()
{
	SolverControl solver_control(1000, 1e-6 * residual_vector.l2_norm());

	SolverGMRES<TrilinosWrappers::MPI::Vector> solver(solver_control);
	TrilinosWrappers::PreconditionSSOR         preconditioner;
	preconditioner.initialize(jacobian_matrix, TrilinosWrappers::PreconditionSSOR::AdditionalData(1.0));

	solver.solve(jacobian_matrix, delta_owned, residual_vector, preconditioner);
	pcout << "   " << solver_control.last_step() << " GMRES iterations" << std::endl;
}

void
MechanicalDisplacement::solve_newton()
{
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
		if (residual_norm > residual_tolerance)
		{
			solve_system();

			solution_owned += delta_owned;
			solution = solution_owned;
		}
		else
		{
			pcout << " < tolerance" << std::endl;
		}

		++n_iter;
	}

	pcout << "===============================================" << std::endl;
}

void
MechanicalDisplacement::output() const
{
	DataOut<dim> data_out;
	data_out.add_data_vector(dof_handler, solution, "u");

	std::vector<unsigned int> partition_int(mesh.n_active_cells());
	GridTools::get_subdomain_association(mesh, partition_int);
	const Vector<double> partitioning(partition_int.begin(), partition_int.end());
	data_out.add_data_vector(partitioning, "partitioning");

	data_out.build_patches();

	const std::string output_file_name = "output-nonlineardiffusion";
	data_out.write_vtu_with_pvtu_record("./", output_file_name, 0, MPI_COMM_WORLD);

	pcout << "Output written to " << output_file_name << "." << std::endl;

	pcout << "===============================================" << std::endl;
}
