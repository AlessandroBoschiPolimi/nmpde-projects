#ifndef MECHANICAL_DISPLACEMENT_HPP
#define MECHANICAL_DISPLACEMENT_HPP

#include <memory>

#include <deal.II/lac/vector.h> // vector
#include <deal.II/base/quadrature.h> // quadrature

#include <deal.II/fe/fe_system.h> // fesystem
#include <deal.II/grid/tria.h> // traingulation

#include <deal.II/dofs/dof_handler.h> // dof handler

#include <deal.II/lac/sparse_matrix.h> // sparse matrix
#include <deal.II/lac/sparsity_pattern.h> // sparsity pattern

// TODO: understand why must use = 0 after virtual function

namespace pde {

using namespace dealii;

class MechanicalDisplacement {

public: static const unsigned int dim = 3;

protected:

    const std::string mesh_file_name;
    const unsigned int r;

    // Neumann conditions
    const std::function<Point<dim>(const Point<dim> &)> h;
    
    Triangulation<dim> mesh;
    const unsigned int num_cells;

    // Finite Element System (USE FE_Q)
    std::unique_ptr<FESystem<dim>> fe;
    
    // Quadratures
    std::unique_ptr<Quadrature<dim>> quadrature;
    std::unique_ptr<Quadrature<dim-1>> quadrature_boundary;

    DoFHandler<dim> dof_handler;

    // Matrices
    SparsityPattern sparsity_pattern;
    SparseMatrix<double> jacobian_matrix;
    Vector<double> residual_vector;
    Vector<double> solution;
    Vector<double> delta;

    virtual void assemble_system() = 0;
    virtual void solve_system() = 0;

public:

    MechanicalDisplacement(
        const std::string mesh_file_name_,
        const unsigned int &r_,
        const std::function<Point<dim>(const Point<dim> &)> &h_,
	const unsigned int num_cells_
    ):
    mesh_file_name(mesh_file_name_),
    r(r_),
    h(h_),
    num_cells(num_cells_)
    {}
    
    virtual void setup() = 0;
    virtual void solve() = 0;
    virtual void output() const = 0;
};

}

#endif // MECHANICAL_DISPLACEMENT_HPP
