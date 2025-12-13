#ifndef TEST_FUNCTIONS_HPP
#define TEST_FUNCTIONS_HPP

#include "MechanicalDisplacement.hpp"
#include <functional>

namespace pde {

/**
* This is a struct that is used to retrieve functions to be applied
* on boundaries as neumann conditions
*/
struct TestNeumannConditions {
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


// TODO: maybe change this so that forcing functions can be dynamically switched

namespace TestForcingFunctions {

static const ForcingTermType null_forcing_term = [](
	const FEValues<dim>& ,
	const FEValuesExtractors::Vector& ,
	const unsigned int ,
	const unsigned int 
) { return 0.0; };

//Forcing term, hardcoded, please dont judge me, otherwise I cannot bend it
// this function f here is f(x,y,z) = (0, 0.02 * x^2, 0) if x > 0
// 				    = (0,0,0)            otherwise 
static const ForcingTermType bend_rod = [](
	const FEValues<dim>& fe_values, 
	const FEValuesExtractors::Vector& displacement, 
	const unsigned int i, 
	const unsigned int q
) {
	return 0.02 * std::pow(fe_values.quadrature_point(q)[0], 2) * 
		fe_values[displacement].value(i, q)[1] * (int)(fe_values.quadrature_point(q)[0] > 0) * fe_values.JxW(q);
};

};

}


#endif // TEST_FUNCTIONS_HPP
