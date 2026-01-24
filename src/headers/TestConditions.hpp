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
			case  2: return 0.2 * (std::sin(p[0] * M_PI * 2));
			default: return 0; // zero displacement along x and y axis
			}
		}
	};
}




namespace TestForcingFunctions {

static const ForcingTermType null_forcing_term = [](const Point<dim> &) {
	Tensor<1, dim> res;
	res[0] = 0;
	res[1] = 0;
	res[2] = 0;	
	return res;
};

// Forcing term, hardcoded, please dont judge me, otherwise I cannot bend it
// this function f here is f(x,y,z) = (0, 0.02 * x^2, 0) if x > 0
// 				    = (0,0,0)            otherwise 
static const ForcingTermType bend_rod = [](const Point<dim> &p) {
	Tensor<1, dim> res;
	res[0] = 0;
	res[1] = 0.02 * std::pow(p[0], 2);
	res[2] = 0;
	return res;
};

static const ForcingTermType cube_squeeze = [](const Point<dim> &p) {
	Tensor<1, dim> res;
	if (p[0] > 0.5) res[0] = -1;
	else            res[0] = +1;
	if (p[1] > 0.5) res[1] = -1;
	else            res[1] = +1;
	if (p[2] > 0.5) res[2] = -1;
	else            res[2] = +1;
	return res;
};
static const ForcingTermType cube_tear = [](const Point<dim> &p) {
	Tensor<1, dim> res;
	if (p[0] > 0.5) res[0] = +1;
	else            res[0] = -1;
	if (p[1] > 0.5) res[1] = +1;
	else            res[1] = -1;
	if (p[2] > 0.5) res[2] = +1;
	else            res[2] = -1;
	return res;
};
static const ForcingTermType cube_squeeze_zm_little = [](const Point<dim> &) {
	Tensor<1, dim> res;
	res[0] = 0;
	res[1] = 0;
	res[2] = -0.2;
	return res;
};
static const ForcingTermType cube_squeeze_zm_lot = [](const Point<dim> &) {
	Tensor<1, dim> res;
	res[0] = 0;
	res[1] = 0;
	res[2] = -2;
	return res;
};

const pde::ForcingTermType choose_forcing_term(std::string func_name);

}

}
