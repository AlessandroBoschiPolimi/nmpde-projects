#ifndef MECHANICAL_DISPLACEMENT_HPP
#define MECHANICAL_DISPLACEMENT_HPP

#include <deal.II/lac/vector.h> // vector
#include <deal.II/base/quadrature.h> // quadrature

#include <deal.II/fe/fe_system.h> // fesystem
#include <deal.II/distributed/fully_distributed_tria.h> // distributed triangulation

#include <deal.II/dofs/dof_handler.h> // dof handler

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h> // sparsity pattern

// ----------- MPI INCLUDES -------------
#include <deal.II/base/conditional_ostream.h>

// ----------- CPP INCLUDES ------------

#include <memory>
#include <unordered_set>

#include "MeshGenerator.hpp"

// TODO: understand why must use = 0 after virtual function

namespace pde {

using namespace dealii;

class MechanicalDisplacement
{
public:
    static const unsigned int dim = 3;

protected:
    const std::unique_ptr<MeshGenerator<dim>> mesh_generator;
    const unsigned int r;

    // Neumann conditions
    const std::function<Point<dim>(const Point<dim> &)> neumann_conds;
    const std::map<types::boundary_id, const Function<dim> *> dirichelet_conds;

    parallel::fullydistributed::Triangulation<dim> mesh;

    // Finite Element System (USE FE_Q)
    std::unique_ptr<FESystem<dim>> fe;
    
    // Quadratures
    std::unique_ptr<Quadrature<dim>> quadrature;
    std::unique_ptr<Quadrature<dim-1>> quadrature_boundary;
    
    // ------------ DOF STUFF ------------------
    DoFHandler<dim> dof_handler;

    // DoFs owned by current process.
    IndexSet locally_owned_dofs;

    // DoFs relevant to the current process (including ghost DoFs).
    IndexSet locally_relevant_dofs;

    // -------- Matrices & Solution Vectors -----------

    // Jacobian matrix.
    TrilinosWrappers::SparseMatrix jacobian_matrix;

    // Residual vector.
    TrilinosWrappers::MPI::Vector residual_vector;

    // Solution increment (without ghost elements).
    TrilinosWrappers::MPI::Vector delta_owned;

    // System solution (without ghost elements).
    TrilinosWrappers::MPI::Vector solution_owned;

    // System solution (including ghost elements).
    TrilinosWrappers::MPI::Vector solution;

    // ------------------ MPI STUFF ---------------------
    // Number of MPI processes.
    const unsigned int mpi_size;
    // This MPI process.
    const unsigned int mpi_rank;
    // Parallel output stream.
    ConditionalOStream pcout;


    std::unordered_set<int> neumann_ids;
    std::string output_filename;

    virtual void assemble_system() = 0;
    virtual void solve_system() = 0;

public:

    MechanicalDisplacement(
            std::unique_ptr<MeshGenerator<dim>> mesh_generator_,
            const unsigned int &r_,
            const std::map<types::boundary_id, const Function<dim> *> boundary_functions_,
            const std::function<Point<dim>(const Point<dim> &)> &neum_funcs_,
            const std::unordered_set<int>& neumann_ids_,
            const std::string& output_filename_
        ) :
        mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
        mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
        pcout(std::cout, mpi_rank == 0),
        mesh_generator(std::move(mesh_generator_)),
        r(r_),
        neumann_conds(neum_funcs_),
        dirichelet_conds(boundary_functions_),
        mesh(MPI_COMM_WORLD),
        neumann_ids(neumann_ids_),
        output_filename(output_filename_)
    {}
    
    virtual void setup() = 0;
    virtual void solve() = 0;
    virtual void output() const = 0;
};

}

#endif // MECHANICAL_DISPLACEMENT_HPP
