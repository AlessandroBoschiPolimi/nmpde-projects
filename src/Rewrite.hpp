#ifndef REWRITE_HPP
#define REWRITE_HPP

#include <deal.II/fe/fe_system.h>

#include <deal.II/dofs/dof_handler.h>

#include <memory>

namespace pde {

using namespace dealii;

class MechanicalDisplacement {

public: static const unsigned int dim = 3;

protected:

    const std::string mesh_file_name;
    const unsigned int r;
    const std::function<Point<dim>(const Point<dim> &)> h;
    
    Triangulation<dim> mesh;

    // Finite Element System (USE FE_Q)
    std::unique_ptr<FESystem<dim>> fe;
    
    // Quadratures
    std::unique_ptr<Quadrature<dim>> quadrature;
    std::unique_ptr<Quadrature<dim-1>> quadrature_boundary;

    DoFHandler<dim> dof_handler;

    // Matrices
    SparsityPattern sparsity_pattern;

    void assemble_system();
    void solve_system();

public:


    MechanicalDisplacement(
        const std::string mesh_file_name_,
        const unsigned int &r_,
        const std::function<Point<dim>(const Point<dim> &)> &h_
    ):
    mesh_file_name(mesh_file_name_),
    r(r_),
    h(h_)
    {}
    
    void setup();
    void solve();
    void output() const;
};

}


#endif // REWRITE_HPP
