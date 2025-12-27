#pragma once
#include "MechanicalDisplacement.hpp"

namespace pde {

class NeoHooke : public MechanicalDisplacement {
protected:
    const double mu;
    const double lambda;
    void assemble_system() override;

public:

    NeoHooke(Config&& config_, const ConditionalOStream pcout_, const unsigned int mpi_rank_, const double mu_, const double lambda_)
        : MechanicalDisplacement(std::move(config_), pcout_, mpi_rank_), mu(mu_), lambda(lambda_)
    {}
};

}
