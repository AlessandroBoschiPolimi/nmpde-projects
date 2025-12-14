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
	const std::string& mesh_file_input_,
        const unsigned int &r_,
	const std::map<types::boundary_id, const Function<dim> *> boundary_functions_,
        const NeumannCondition &neum_conds_,
	const ForcingTermType forcing_term_,
	const double mu_,
    	const double lambda_
    ) :
        MechanicalDisplacement(mesh_file_input_, r_, boundary_functions_, neum_conds_, forcing_term_),
        mu(mu_),
        lambda(lambda_)
    {}

    void setup() override;
    void solve() override;
    void output() const override;
};

}
