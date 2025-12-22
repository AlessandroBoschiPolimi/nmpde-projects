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
constexpr auto bowl_pull_out = [](const Point<dim> &p) {
	const double tau = TestNeumannConditions::parameter;
	double xx = p[0] * p[0], yy = p[1] * p[1], zz = p[2] * p[2];
	return tau / std::sqrt(xx + yy + zz / 81.0) *
		Point<dim>(p[0], p[1], p[2] / 9.0); // vector normal to the outer surface of the acorn;
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
    if (choice.starts_with("bowl"))
    	return functions::bowl_pull_out;
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
