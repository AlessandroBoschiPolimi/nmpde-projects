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


namespace TestForcingFunctions {

static const ForcingTermType null_forcing_term = [](const Point<dim> &) { return 0.0; };

//Forcing term, hardcoded, please dont judge me, otherwise I cannot bend it
// this function f here is f(x,y,z) = (0, 0.02 * x^2, 0) if x > 0
// 				    = (0,0,0)            otherwise 
static const ForcingTermType bend_rod = [](const Point<dim> &p) {
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

constexpr auto rod_pull = [](const Point<dim> &p) {
    const double tau = TestNeumannConditions::parameter;
    return Point<dim>(0, 0, tau * (int)(p[2] < -0.2));
};


constexpr auto cube_push = [](const Point<dim> &p) {
    constexpr double small_tol = 1e-13;
    const double tau = TestNeumannConditions::parameter;

    if (std::abs(p[0]) < small_tol)
	    return Point<dim>(tau, 0, 0);
    else if (std::abs(p[0]) > (1 - small_tol))
	    return Point<dim>(-tau, 0, 0);
    else if (std::abs(p[1]) < small_tol)
	    return Point<dim>(0, tau, 0);
    else if (std::abs(p[1]) > (1 - small_tol))
	    return Point<dim>(0, -tau, 0);
    else if (std::abs(p[2]) < small_tol)
	    return Point<dim>(0, 0, tau);
    else if (std::abs(p[2]) > (1 - small_tol))
	    return Point<dim>(0, 0, -tau);
    else
	    return Point<dim>(0, 0, 0);
};

constexpr auto cube_push_z = [](const Point<dim> &p) {
    if(p[1] > 0.4 &&  p[1] < 0.6 && p[0] > 0.4 && p[0] < 0.6)
	    return Point<dim>(0, 0, TestNeumannConditions::parameter);
    else 
    	return Point<dim>(0,0,0);	
};

constexpr auto force_cube_push_z = [](const Point<dim> &p) {
    return Point<dim>(0,0,TestNeumannConditions::parameter);
};

// Applies a force in the normal direction w.r.t. surface
// of an ellipsoid with semi-axes with the ratio 1:1:3. OUTGOING
constexpr auto bowl_pull_out = [](const Point<dim> &p) {
	const double tau = TestNeumannConditions::parameter;
	double xx = p[0] * p[0], yy = p[1] * p[1], zz = p[2] * p[2];
    if (xx + yy + zz / 9 > 0.9)
	return tau / std::sqrt(xx + yy + zz / 81.0) *
		Point<dim>(p[0], p[1], p[2] / 9.0); // vector normal to the outer surface of the acorn;
	return Point<dim>(0, 0, 0);
};

}

}

#endif // TEST_FUNCTIONS_HPP
