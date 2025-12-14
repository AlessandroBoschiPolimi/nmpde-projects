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
        std::unique_ptr<MeshGenerator<dim>> mesh_generator_,
        const unsigned int &r_,
	const std::map<types::boundary_id, const Function<dim> *> boundary_functions_,
        const std::function<Point<dim>(const Point<dim> &)> &neum_funcs_,
    	const double mu_,
    	const double lambda_,
        const std::unordered_set<int>& neumann_ids_,
        const std::string& output_filename_
    ) :
        MechanicalDisplacement(std::move(mesh_generator_), r_, boundary_functions_, neum_funcs_, neumann_ids_, output_filename_),
        mu(mu_),
        lambda(lambda_)
    {}

    void setup() override;
    void solve() override;
    void output() const override;
};

}
