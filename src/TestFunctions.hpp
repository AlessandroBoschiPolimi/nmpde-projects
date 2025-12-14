#ifndef TEST_FUNCTIONS_HPP
#define TEST_FUNCTIONS_HPP

#include "MechanicalDisplacement.hpp"
#include <functional>

namespace pde {

struct TestFunctions {
    // ----- PARAMETERS -------
    static constexpr unsigned int dim = MechanicalDisplacement::dim;
    // ------- VARIABLES ----

    // Copying from tau now but might change it in the future
    static double parameter;

    // ------- FUNCTIONS ------

    static void initialize(double tau_0 = 0.5);
    static const std::function<Point<dim> (const Point<dim> &)> 
	choose_neumann_function(std::string choice);
    };
}


#endif // TEST_FUNCTIONS_HPP
