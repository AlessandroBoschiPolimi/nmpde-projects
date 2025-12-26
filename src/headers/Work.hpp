#include <optional>
#include <vector>
#include <string>
#include <unordered_set>
#include <variant>
#include <filesystem>

#include "defs.hpp"

namespace pde {

struct Work {
    enum MaterialType { NeoHooke, Guccione };
    MaterialType material;
    enum GeometryType { File, Cube, Rod };
    GeometryType geometry;

    // GeometryType::File: file path
    // GeometryType::Cube: mesh refinement
    std::optional<std::variant<std::filesystem::path, unsigned int>> mesh_param;

    std::string output_filename;
    int iterations = 0;

    bool newton_damping = false;
    double newton_scaling = 1.0;

    std::unordered_set<int> N_values;
    std::string N_label;
    std::string N_data;

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
};

std::vector<Work> parse_file(const std::filesystem::path& path);
}
