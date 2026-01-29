#pragma once

#include <deal.II/base/conditional_ostream.h>

#include <optional>
#include <vector>
#include <string>
#include <unordered_set>
#include <variant>
#include <filesystem>

#include "defs.hpp"

namespace pde {

struct Work {
    std::string name;
    bool passed;

    enum MaterialType { NeoHooke, Guccione };
    MaterialType material;
    enum GeometryType { File, Cube, Rod };
    GeometryType geometry;

    // GeometryType::File: file path
    // GeometryType::Cube: mesh refinement
    std::optional<std::variant<std::filesystem::path, unsigned int>> mesh_param;

    std::string output_filename;
    // GMRES maximum iterations per newton step
    int iterations = 0;

    double newton_damping = 1.0;

    // name of the forcing term function to use
    std::string forcing_term;

    std::unordered_set<int> neumann_ids; // boundary ids of the neumann condition
    std::string neumann_func;            // name of the neumann function
    std::string neumann_func_param;      // parameter of the neumann function

    struct DEntry {
        int value;
        std::string function;
    };
    std::vector<DEntry> dirichlet_conditions;

    struct NeoHookeData
    {
        double mu, lambda;
    };
    struct GuccioneData
    {
        double c;
        std::array<double, 9> b;
        std::array<Point<dim>, dim> aniso_fun_points;
    };
    std::variant<NeoHookeData, GuccioneData> problem_params;

    bool time = false;
};

std::vector<Work> parse_file(const std::filesystem::path& path, const ConditionalOStream &pcout);
}
