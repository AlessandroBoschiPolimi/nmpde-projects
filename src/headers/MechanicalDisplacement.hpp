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
#include <filesystem>

#include "MeshGenerator.hpp"
#include "defs.hpp"

namespace pde {

class MechanicalDisplacement
{

protected:
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
    const unsigned int mpi_rank;
    const ConditionalOStream pcout;


    virtual void assemble_system() = 0;
    virtual void solve_system() = 0;

public:
    struct Config
    {
        const int iterations = 10000;
        const std::string output_filename;

        std::unique_ptr<MeshGenerator<dim>> mesh_generator;
        const unsigned int r;

        const bool newton_damping = false;
        const double newton_scaling = 0.2;

        // ---------- NEUMANN CONDITIONS ----------
        const std::function<Point<dim>(const Point<dim> &)> neumann_conds;
        const std::unordered_set<int> neumann_ids;
        
        // ---------- DIRICHELET CONDITIONS -----------
        const std::map<types::boundary_id, const Function<dim> *> dirichelet_conds;
        
        // --------- FORCING TERM -------------
        const ForcingTermType forcing_term;

        Config(const int iterations_, const std::string& output_filename_, 
                std::unique_ptr<MeshGenerator<dim>>&& mesh_generator_, const unsigned int r_,
                const bool newton_damping_, const double newton_scaling_,
                const std::function<Point<dim>(const Point<dim> &)>& neumann_conds_,
                const std::unordered_set<int>& neumann_ids_, 
                const std::map<types::boundary_id, const Function<dim> *>& dirichelet_conds_,
                const ForcingTermType& forcing_term_) :
            iterations(iterations_), output_filename(output_filename_),
            mesh_generator(std::move(mesh_generator_)), r(r_),
            newton_damping(newton_damping_), newton_scaling(newton_scaling_),
            neumann_conds(neumann_conds_), neumann_ids(neumann_ids_),
            dirichelet_conds(dirichelet_conds_), forcing_term(forcing_term_)
        {}

        Config(Config&& other) :
            iterations(other.iterations), output_filename(other.output_filename),
            mesh_generator(std::move(other.mesh_generator)), r(other.r),
            newton_damping(other.newton_damping), newton_scaling(other.newton_scaling),
            neumann_conds(other.neumann_conds), neumann_ids(other.neumann_ids),
            dirichelet_conds(other.dirichelet_conds), forcing_term(other.forcing_term)
        {}
    };

protected:
    Config config;

public:

    MechanicalDisplacement(
            Config&& config_,
	        const ConditionalOStream pcout_,
	        const unsigned int mpi_rank_
    ) :
	    mesh(MPI_COMM_WORLD),
        mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
	    mpi_rank(mpi_rank_),
	    pcout(pcout_),
        config(std::move(config_))
    {}
    
    virtual void setup() = 0;
    virtual void solve() = 0;
    virtual void output() const = 0;
};

}

#endif // MECHANICAL_DISPLACEMENT_HPP
