#include "MechanicalDisplacement.hpp"

namespace pde {

class NeoHooke : public MechanicalDisplacement {
protected:
    const double mu;
    const double lambda;
    void assemble_system() override;
    void solve_system() override;

    void apply_dirchlet_to_initial_solution();
    void apply_zero_dirchlet_to_newton_update();

public:

    NeoHooke(Config&& config_, const ConditionalOStream pcout_, const unsigned int mpi_rank_, const double mu_, const double lambda_)
        : MechanicalDisplacement(std::move(config_), pcout_, mpi_rank_), mu(mu_), lambda(lambda_)
    {}

    void setup() override;
    void solve() override;
    void output() const override;
};

}
