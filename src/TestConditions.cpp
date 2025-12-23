#include "headers/TestConditions.hpp"
#include <cmath>

namespace pde {

// ----- variables declarations  ------
// without this it would give undefined reference

double TestNeumannConditions::parameter;

// ---------------- FUNCTIONS --------------------------

namespace functions {

/// ---------------- CUBE SECTION -----------------

/**
* This is the function that we used until now to
* pull the cube on x and y
*/
const auto cube_pull = [&tau = TestNeumannConditions::parameter](const Point<dim> &p) {
    constexpr double small_tol = 1e-13;

    if (std::abs(p[0]) < small_tol)
	    return Point<dim>(-tau, 0, 0);
    else if (std::abs(p[0]) > (1 - small_tol))
	    return Point<dim>(tau, 0, 0);
    else if (std::abs(p[1]) < small_tol)
	    return Point<dim>(0, -tau, 0);
    else if (std::abs(p[1]) > (1 - small_tol))
	    return Point<dim>(0, tau, 0);
    else if (std::abs(p[2]) < small_tol)
	    return Point<dim>(0, 0, -tau);
    else if (std::abs(p[2]) > (1 - small_tol))
	    return Point<dim>(0, 0, tau);
    else
	    return Point<dim>(0, 0, 0);

};

const auto cube_push = [&tau = TestNeumannConditions::parameter](const Point<dim> &p) {
    constexpr double small_tol = 1e-13;

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


/// --------------- BOWL SECTION ----------------

// Applies a force with modulus 'tau' in the normal direction w.r.t. surface
// of an ellipsoid with semi-axes with the ratio 1:1:3. OUTGOING
// This is the outer surface of our acorn
constexpr auto bowl_pull_out = [](const Point<dim> &p) {
	const double tau = TestNeumannConditions::parameter;
	double xx = p[0] * p[0], yy = p[1] * p[1], zz = p[2] * p[2];
	return tau / std::sqrt(xx + yy + zz / 81.0) *
		Point<dim>(p[0], p[1], p[2] / 9.0); // vector normal to the outer surface of the acorn;
};

// Applies a force with modulus 'tau' in the normal direction w.r.t. surface
// of an ellipsoid with semi-axes with the ratio 7:7:17. OUTGOING
// This is the force used in paper 5 problem 2, inside the acorn.
constexpr auto bowl_push_in = [](const Point<dim> &p) {
	const double tau = TestNeumannConditions::parameter;
	double xx = p[0] * p[0], yy = p[1] * p[1], zz = p[2] * p[2];
    constexpr double a2 = 7.0 * 7.0, c2 = 17.0 * 17.0;
    constexpr double a4 = a2 * a2, c4 = c2 * c2;
	return tau / std::sqrt(xx / a4 + yy / a4 + zz / c4) * Point<dim>(p[0] / a2, p[1] / a2, p[2] / c2);
};

/// -------------  ROD SECTION ----------------

constexpr auto rod_pull = [](const Point<dim> &p) {
    const double tau = TestNeumannConditions::parameter;
    return Point<dim>(0, 0, tau);
};



}


const std::function<Point<dim> (const Point<dim> &)> 
	TestNeumannConditions::choose_neumann_function(std::string choice)
{
    // Here I simply started defining different functions
    // The more models and boundary conditions we apply the better it is.
    if (choice == "bowl")
    	return functions::bowl_pull_out;
    if (choice == "bowl_ref")
    	return functions::bowl_push_in;
    if (choice == "cube_pull")
	    return functions::cube_pull;
    if (choice == "rod_pull")
    	return functions::rod_pull;
    if (choice == "cube_push")
    	return functions::cube_push;

    throw std::runtime_error("Unknown Neumann Condition: " + choice);
}


void TestNeumannConditions::initialize(double tau_0) {
    parameter = tau_0;
}

}
