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
    int iterations = 0;

    double newton_damping = 1.0;

    std::string forcing_term;

    std::unordered_set<int> N_values; // boundary ids of the neumann condition
    std::string N_label; // name of the neumann function
    std::string N_data;  // parameter of the neumann function

    struct DEntry {
        int value;
        std::string function;
    };
    std::vector<DEntry> D_entries;

    struct NeoHookeData
    {
        double C, lambda;
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

std::vector<Work> parse_file(const std::filesystem::path& path,
			     const ConditionalOStream &pcout);
}
