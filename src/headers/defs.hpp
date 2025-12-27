#pragma once

#include <deal.II/base/point.h>

namespace pde {

using namespace dealii;
static constexpr unsigned int dim = 3;
using ForcingTermType = std::function<double (const Point<dim>& p)>;
using AnisotropicFunctionType = std::function<std::array<Point<dim>, dim> (const Point<dim> &)>;

}