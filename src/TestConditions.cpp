#include "headers/TestConditions.hpp"
#include <cmath>

namespace pde {

// ----- variables declarations  ------
// without this it would give undefined reference

double TestNeumannConditions::parameter;

// TODO: Implement this function

const std::function<Point<dim> (const Point<dim> &)> 
	TestNeumannConditions::choose_neumann_function(std::string choice)
{
    // Here I simply started defining different functions
    // The more models and boundary conditions we apply the better it is.
    if (choice.starts_with("bowl")) {
	return functions::bowl_pull_out;
    }
    if(choice == "cube_pull") {
	return functions::cube_pull;
    }
    if(choice == "cube_push") {
	return functions::cube_push;
    }
    if(choice == "cube_push_x") {
		std::cout << "Running cube_push_x" << std::endl;
	return functions::cube_push_x;
    }
    throw std::runtime_error("Unknown Neumann Condition: " + choice);
}


void TestNeumannConditions::initialize(double tau_0) {
    parameter = tau_0;
}

}
