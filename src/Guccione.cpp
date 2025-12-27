#include "headers/Guccione.hpp"

#define GAMBA_DEBUG false

using namespace pde;


void Guccione::assemble_system() {
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
		update_values | update_quadrature_points | update_normal_vectors |
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


    for (const auto &cell : dof_handler.active_cell_iterators()) {
		if (!cell->is_locally_owned())
			continue;

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
            // evaluate anisotropic function
            const std::array<Point<dim>,dim> fns = aniso_fun(fe_values.quadrature_point(q));
	    	//Get a unit matrix, cause we will need it
            Tensor<2,dim> I({{1,0,0},{0,1,0},{0,0,1}});
			//Build the deformation gradient
			Tensor<2,dim> F = I + solution_gradient_loc[q];
            // build E
            Tensor<2, dim> E = 0.5 * (transpose(F)*(F) - I);
            // build D and B
            Tensor<4, dim> D;
            Tensor<2, dim> B;
            for (unsigned int aniso_vect_1 = 0; aniso_vect_1 < dim; ++aniso_vect_1) {
                auto m = fns[aniso_vect_1];
                for (unsigned int aniso_vect_2 = 0; aniso_vect_2 < dim; ++aniso_vect_2) {
                    auto n = fns[aniso_vect_2];
                    const double b = param_b[aniso_vect_1 * dim + aniso_vect_2];
                    for (unsigned int i = 0; i< dim; i++) {
                    	for (unsigned int j = 0; j< dim; j++) {
							if (i != j)
								B[i][j] +=  2 * b * ((E * m) * n) * (m[i] * n[j] + m[j] * n[i]);
							else
								B[i][j] +=  2 * b * ((E * m) * n) * (m[i] * n[j]);
							for (unsigned int k = 0; k< dim; k++) {
								for (unsigned int l = 0; l< dim; l++) {
									if ((k != l) && (i != j))
										D[k][l][i][j] +=  2 * b * (m[i] * n[j] + m[j] * n[i])
											* (m[k] * n[l] + m[l] * n[k]);
									else if ((k == l) && (i != j))
										D[k][l][i][j] += 2 * b * (m[i] * n[j] + m[j] * n[i])
													* (m[k] * n[l]);
									else if ((k != l) && (i == j))
										D[k][l][i][j] +=  2 * b * (m[i] * n[j])
													* (m[k] * n[l] + m[l] * n[k]);
									else if ((k == l) && (i == j))
										D[k][l][i][j] +=  2 * b * (m[i] * n[j])
													* (m[k] * n[l]);
								}
							}
						}
                    }
                }
            }
			// compute Q
			double Q = 0;
            for (unsigned int aniso_vect_1 = 0; aniso_vect_1 < dim; ++aniso_vect_1) {
                auto m = fns[aniso_vect_1];
                for (unsigned int aniso_vect_2 = 0; aniso_vect_2 < dim; ++aniso_vect_2) {
                    auto n = fns[aniso_vect_2];
                    const double b = param_b[aniso_vect_1 * dim + aniso_vect_2];
					Q += b * std::pow((E * m) * n, 2);
				}
	    	}
            // build P
	    	Tensor<2, dim> P;
            Tensor<2, dim> Ft = transpose(invert(F));
	    	double det = std::abs(determinant(F));
	    	double B_modulus = 0;
	    	P = param_c / 2.0 * std::exp(Q) * F * B + B_modulus / 2.0 * (det * std::log(det) + det - 1) * Ft;
            // compute dP/dF
			Tensor<4, dim> dPdF;
			for (unsigned int k = 0; k < dim; ++k){
				for (unsigned int l = 0; l < dim; ++l){
					for (unsigned int i = 0; i < dim; ++i){
						for (unsigned int j = 0; j < dim; ++j){
							dPdF[k][l][i][j] -= B_modulus / 2.0 * Ft[i][l] * Ft[k][j]  * (det * std::log(det) + det - 1);
							dPdF[k][l][i][j] += B_modulus / 2.0 * Ft[k][l] * Ft[i][j] * det * (std::log(det) + 2);
							if (k == i){
								dPdF[k][l][i][j] += param_c / 2.0 * std::exp(Q) *
									(B[l][j]);
							}
							dPdF[k][l][i][j] += param_c / 2.0 * std::exp(Q) *
										(F * B)[k][l] * (F * B)[i][j];
							for (unsigned int n = 0; n < dim; n++){
								for (unsigned int a = 0; a < dim; a++){
									dPdF[k][l][i][j] += param_c / 2.0 * std::exp(Q) * F[k][a] * F[i][n] * (D[a][l][n][j]);
								}
							}
						}
					}
				}
			}
			for ( const unsigned int i : fe_values.dof_indices() ) {
				// base to describe v
				const Tensor<2,dim> phi_i_grad = fe_values[displacement].gradient(i, q);
				for ( const unsigned j : fe_values.dof_indices() ) {
					// base to describe delta
					const Tensor<2, dim> phi_j_grad = fe_values[displacement].gradient(j,q);
					//TODO:Check the signs
					//(dPdF:grad(delta)):grad(v)
					//TODO: Does this contract the last two indeces? Is the inermediate right?
					Tensor<2, dim> intermediate;
					for (unsigned int a = 0; a < dim; a++){
						for (unsigned int b = 0; b < dim; b++){
							for (unsigned int c = 0; c < dim; c++){
								for (unsigned int d = 0; d < dim; d++){
									intermediate[a][b] += dPdF[a][b][c][d] * phi_j_grad[c][d];
								}
							}
						}
					}
					cell_matrix(i, j) += double_contract<0,0,1,1>(intermediate, phi_i_grad) * fe_values.JxW(q) ;
				}

				cell_rhs(i) -= double_contract<0,0,1,1>(P, phi_i_grad) * fe_values.JxW(q);
			}
		}
		
		if (cell->at_boundary()) {
			for (unsigned int face_number = 0; 
				face_number < cell->n_faces();
				face_number++) {

				const unsigned int bound_id = cell->face(face_number)->boundary_id();

				if (!cell->face(face_number)->at_boundary())
					continue;
				if (!config.neumann_ids.contains(bound_id))
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

						// cell_rhs(i) +=
						// 	-10000 *
						// 	fe_values_boundary.normal_vector(q) *
						// 	fe_values_boundary[displacement].value(i, q) * 
						// 	fe_values_boundary.JxW(q);

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
    
    #if GAMBA_DEBUG
	pcout << "system_rhs " << residual_vector << std::endl;
    #endif


    {
		std::map<types::global_dof_index, double> boundary_values;

		VectorTools::interpolate_boundary_values(dof_handler,
					config.dirichelet_conds,
					boundary_values);

		MatrixTools::apply_boundary_values(boundary_values, jacobian_matrix, delta_owned, residual_vector, false);
    }
}

