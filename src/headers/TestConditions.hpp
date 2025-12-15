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
	const Point<dim> &
) { return 0.0; };

//Forcing term, hardcoded, please dont judge me, otherwise I cannot bend it
// this function f here is f(x,y,z) = (0, 0.02 * x^2, 0) if x > 0
// 				    = (0,0,0)            otherwise 
static const ForcingTermType bend_rod = [](
	const Point<dim> &p 
) {
	return 0.02 * std::pow(p[0], 2);
};

}

// ---------------- FUNCTIONS --------------------------

namespace functions {

// TODO: add explanation of functions

/**
* This is the function that we used until now to
* pull the cube on x and y
*/
constexpr auto cube_pull = [](const Point<dim> &p) {
    constexpr double small_tol = 1e-13;
    const double tau = TestNeumannConditions::parameter;

    if (std::abs(p[0]) < small_tol && std::abs(p[1]) > small_tol)
	    return Point<dim>(-tau, 0, 0);
    else if (std::abs(p[0]) > (1 - small_tol) && std::abs(p[1]) > small_tol)
	    return Point<dim>(tau, 0, 0);
    else if (std::abs(p[1]) < small_tol && std::abs(p[0]) > small_tol)
	    return Point<dim>(0, -tau, 0);
    else if (std::abs(p[1]) > (1 - small_tol) && std::abs(p[0]) > small_tol)
	    return Point<dim>(0, tau, 0);
    else
	    return Point<dim>(0, 0, -tau);
};

// FIXME: This function doesn't work
// the simulation doesn't converge
/**
* This function should bend the rod as a horseshoe
* it applies a force along z based on the x coordinate 
* and proportional to tau
*/
constexpr auto rod_bend = [](const Point<dim> &p) {
    const double tau = TestNeumannConditions::parameter;
    // TODO: check that this function is in fact correct 
    // and it's not the cause of the divergence
    return Point<dim>(0, 0, tau * p[0]);
};

constexpr auto bowl_pull_out = [](const Point<dim> &p) {
	static const auto aa = 1 * 1, bb = 1 * 1, cc = 3 * 3;
	const double tau = TestNeumannConditions::parameter;
	auto xx = p[0] * p[0], yy = p[1] * p[1], zz = p[2] * p[2];
	if (xx / aa + yy / bb + zz / cc > 0.9)
	    return tau / std::sqrt(64 * xx + 64 * yy + 64 / 81.0 * zz) *
		Point<dim>(8 * p[0], 8 * p[1], 8 / 9.0 * p[2]); // vector normal to the outer surface of the acorn;
	return Point<dim>(0, 0, 0);
};

}

}

#endif // TEST_FUNCTIONS_HPP
