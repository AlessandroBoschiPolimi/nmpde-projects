#include <optional>
#include <vector>
#include <string>
#include <unordered_set>
#include <variant>

#include "defs.hpp"

namespace pde {

struct Work {
    enum MaterialType { NeoHooke, Guccione };
    MaterialType material;
    enum GeometryType { File, Cube, Rod };
    GeometryType geometry;
    std::optional<std::string> filename; // only for GeometryType::File

    std::string output_filename;
    int iterations = 0;

    std::unordered_set<int> N_values;
    std::string N_label;
    std::string N_data;

    struct DEntry {
        int value;
        // expand with different function on D boundary
    };
    std::vector<DEntry> D_entries;

    struct NeoHookeData
    {
        double C, lambda;
    };
    struct GuccioneData
    {
        double c, alpha;
        std::array<double, 9> b;
        std::array<Point<dim>, dim> aniso_fun_points;
    };
    std::variant<NeoHookeData, GuccioneData> problem_params;

    /*
    double lambda_param;
    double C_param, alpha_param;
    std::array<double, 9> B_param;
    std::array<Point<dim>, dim> aniso_fun_points;
    */
};

std::vector<Work> parse_file(const std::string& path);
}
