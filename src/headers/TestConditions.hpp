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
	static const std::function<Point<dim>(const Point<dim> &)> choose_neumann_function(std::string func_name);
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
			case  2: return 0.05 * (std::sin(p[0] * M_PI * 3));
			default: return 0; // zero displacement along x and y axis
			}
		}
	};
}




namespace TestForcingFunctions {

const ForcingTermType choose_forcing_term(std::string func_name);

}

}
