#include "MechanicalDisplacement.hpp"


namespace pde {

class NeoHooke : public MechanicalDisplacement {
protected:
    const double mu;
    const double lambda;
    void assemble_system() override;
    void solve_system() override;

public:

    NeoHooke(
        const std::string mesh_file_name_,
        const unsigned int &r_,
	const std::map<types::boundary_id, const Function<dim> *> boundary_functions_,
        const std::function<Point<dim>(const Point<dim> &)> &neum_funcs_,
	const unsigned int num_cells_,
	const double mu_,
	const double lambda_
    ):
    MechanicalDisplacement(mesh_file_name_, r_, boundary_functions_, neum_funcs_,  num_cells_),
    mu(mu_),
    lambda(lambda_)
    {}

    void setup() override;
    void solve() override;
    void output() const override;
};

}
