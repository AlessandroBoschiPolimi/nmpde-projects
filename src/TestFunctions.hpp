#ifndef TEST_FUNCTIONS_HPP
#define TEST_FUNCTIONS_HPP

#include "MechanicalDisplacement.hpp"
#include <functional>

namespace pde {

/**
* This is a struct that is used to retrieve functions to be applied
* on boundaries as neumann conditions
*/
struct TestFunctions {
    // ----- PARAMETERS -------
    static constexpr unsigned int dim = MechanicalDisplacement::dim;
    // ------- VARIABLES ----

    // Copying from tau now but might change it in the future
    static double parameter;

    // ------- FUNCTIONS ------

    // Initializes parameters
    static void initialize(double tau_0 = 0.5);
    // Retrieves the function based on the name
    static const std::function<Point<dim> (const Point<dim> &)> 
	choose_neumann_function(std::string func_name);
    };
}


#endif // TEST_FUNCTIONS_HPP
