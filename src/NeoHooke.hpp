#include "MechanicalDisplacement.hpp"


namespace pde {

class NeoHooke : public MechanicalDisplacement {
protected:
    void assemble_system() override;
    void solve_system() override;

public:

    NeoHooke(
        const std::string mesh_file_name_,
        const unsigned int &r_,
        const std::function<Point<dim>(const Point<dim> &)> &h_,
	const unsigned int num_cells_
    ):
    MechanicalDisplacement(mesh_file_name_, r_, h_, num_cells_)
    {}

    void setup() override;
    void solve() override;
    void output() const override;
};

}
