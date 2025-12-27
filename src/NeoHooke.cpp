#include "headers/NeoHooke.hpp"

#define GAMBA_DEBUG false

using namespace pde;


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
			if (determinant_displacement <= 0)
  				throw std::runtime_error("det(F) <= 0");

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
