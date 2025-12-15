#ifndef GUCCIONE_HPP
#define GUCCIONE_HPP

#include "MechanicalDisplacement.hpp"

namespace pde {


class Guccione : public MechanicalDisplacement {
protected:
    const double param_c, param_alpha;
    const std::array<double, 9> param_b;
    // REturns the f, n, s vectors
    const std::function<std::array<Point<dim>, dim>(const Point<dim> &)> &aniso_fun; 

    void assemble_system() override;
    void solve_system() override;

public:

    Guccione(
        std::unique_ptr<MeshGenerator<dim>> mesh_generator_,
        const unsigned int &r_,
	    const std::map<types::boundary_id, const Function<dim> *> boundary_functions_,
        const std::function<Point<dim>(const Point<dim> &)> &neum_funcs_,
        const std::unordered_set<int>& neumann_ids_,
	    const ForcingTermType forcing_term_,
        const std::string& output_filename_,
        const int iterations_,
	    const ConditionalOStream pcout_,
	    const unsigned int mpi_rank_,
	    const double param_c_,
        const std::array<double, 9> param_b_,
        const AnisotropicFunctionType &aniso_fun_,
        const double param_alpha_
    ) :
        MechanicalDisplacement(
	    std::move(mesh_generator_), r_, 
	    boundary_functions_, neum_funcs_, neumann_ids_,
	    forcing_term_, output_filename_, iterations_, pcout_, mpi_rank_
	),
        param_c(param_c_),
        param_b(param_b_),
        aniso_fun(aniso_fun_),
        param_alpha(param_alpha_)
    {}

    void setup() override;
    void solve() override;
    void output() const override;
};

}

#endif
