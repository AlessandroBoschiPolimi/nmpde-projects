#include "TestConditions.hpp"

#define GAMBA_DEBUG true

namespace pde {

// ----- variables declarations  ------
// without this it would give undefined reference

double TestNeumannConditions::parameter;

namespace functions {

// TODO: add explanation of functions

/**
* This is the function that we used until now to
* pull the cube on x and y
*/
constexpr auto cube_pull = [](const Point<dim> &p) {
    constexpr double small_tol = 1e-13;
    const double tau = TestNeumannConditions::parameter;

    #if GAMBA_DEBUG
    #endif
    if (std::abs(p[0]) < small_tol && std::abs(p[1]) > small_tol)
	    return Point<dim>(-tau, 0, 0);
    else if (std::abs(p[0]) > (1 - small_tol) && std::abs(p[1]) > small_tol)
	    return Point<dim>(tau, 0, 0);
    else if (std::abs(p[1]) < small_tol && std::abs(p[0]) > small_tol)
	    return Point<dim>(0, -tau, 0);
    else if (std::abs(p[1]) > (1 - small_tol) && std::abs(p[0]) > small_tol)
	    return Point<dim>(0, tau, 0);
    else
	    return Point<dim>(0, 0, 0);
};

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

}

// TODO: Implement this function

const std::function<Point<dim> (const Point<dim> &)> 
	TestNeumannConditions::choose_neumann_function(std::string choice)
{
    // Here I simply started defining different functions
    // The more models and boundary conditions we apply the better it is.
    if(choice.starts_with("cube")) {
	return functions::cube_pull;
    } else {
	#if GAMBA_DEBUG
	std::cout << "Calling rod bend" << std::endl;
	#endif
	return functions::rod_bend;
    }
}


void TestNeumannConditions::initialize(double tau_0) {
    parameter = tau_0;
}

}
