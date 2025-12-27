#pragma once

#include <deal.II/base/point.h>
#include <deal.II/base/function.h>

#include <functional>
#include "defs.hpp"

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

namespace TestDirichletConditions {
    class SinXYFunction : public Function<dim>
    {
    public:
        SinXYFunction() : Function<dim>(dim) {};
    
        virtual double value(const Point<dim> &p, const unsigned int component = 0) const override
        { 
            switch (component)
            {
            case  2: return 0.2 * (std::sin(p[0]));
            default: return 0; // zero displacement along x and y axis
            }
        }
    };
}




namespace TestForcingFunctions {

static const ForcingTermType null_forcing_term = [](const Point<dim> &) { return 0.0; };

// Forcing term, hardcoded, please dont judge me, otherwise I cannot bend it
// this function f here is f(x,y,z) = (0, 0.02 * x^2, 0) if x > 0
// 				    = (0,0,0)            otherwise 
static const ForcingTermType bend_rod = [](const Point<dim> &p) {
	return 0.02 * std::pow(p[0], 2);
};

}

}
