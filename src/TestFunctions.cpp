#include <TestFunctions.hpp>

#define GAMBA_DEBUG false

namespace pde {

// ----- variables declarations  ------
// without this it would give undefined reference

double TestFunctions::parameter;

constexpr unsigned int dim = TestFunctions::dim;

namespace functions {

// TODO: add explanation of functions

constexpr auto cube_pull = [](const Point<dim> &p) {
    constexpr double small_tol = 1e-13;
    const double tau = TestFunctions::parameter;

    #if GAMBA_DEBUG
    std::cout << "parameter " << tau << "boh " << TestFunctions::parameter << std::endl;
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

}

// TODO: Implement this function

const std::function<Point<dim> (const Point<dim> &)> 
	TestFunctions::choose_neumann_function(std::string)
{
	return functions::cube_pull;
}


void TestFunctions::initialize(double tau_0) {
    parameter = tau_0;
}

}
