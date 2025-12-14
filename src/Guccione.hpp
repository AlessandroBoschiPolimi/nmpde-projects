#include "MechanicalDisplacement.hpp"

namespace pde {

class Guccione : public MechanicalDisplacement {
protected:
    const double param_c;
    const std::array<double, 9> param_b;
    // REturns the f, n, s vectors
    const std::function<std::array<Point<dim>, dim>(const Point<dim> &)> &aniso_fun; 

    void assemble_system() override;
    void solve_system() override;

public:

    Guccione(
	const std::string mesh_file_input_,
        const unsigned int &r_,
	const std::map<types::boundary_id, const Function<dim> *> boundary_functions_,
	const NeumannCondition &neum_conds_,
	const ForcingTermType forcing_term_,
	const double param_c_,
        const std::array<double, 9> param_b_,
        const std::function<std::array<Point<dim>, dim>(const Point<dim> &)> &aniso_fun_):
        MechanicalDisplacement(mesh_file_input_, r_, boundary_functions_, neum_conds_, forcing_term_),
        param_c(param_c_),
        param_b(param_b_),
        aniso_fun(aniso_fun_)
    {}

    void setup() override;
    void solve() override;
    void output() const override;
};

}
