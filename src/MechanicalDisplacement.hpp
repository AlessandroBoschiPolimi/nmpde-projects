#ifndef REWRITE_HPP
#define REWRITE_HPP

#include <deal.II/fe/fe_system.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/sparse_matrix.h>
#include <memory>

// TODO: understand why must use = 0 after virtual function

namespace pde {

using namespace dealii;

class MechanicalDisplacement {

public: static const unsigned int dim = 3;

protected:

    const std::string mesh_file_name;
    const unsigned int r;
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
    SparseMatrix<double> system_matrix;
    Vector<double> system_rhs;
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


#endif // REWRITE_HPP
