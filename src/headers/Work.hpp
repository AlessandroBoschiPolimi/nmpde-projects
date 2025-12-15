#include <optional>
#include <vector>
#include <string>
#include <unordered_set>

#include "defs.hpp"

namespace pde {

struct Work {
    enum MaterialType { NeoHooke, Guccione };
    MaterialType material;
    enum GeometryType { File, Cube, Rod };
    GeometryType geometry;
    std::optional<std::string> filename; // only for file
    std::string output_filename;
    int iterations = 0;

    std::unordered_set<int> N_values;
    std::string N_label;
    std::string N_data;

    double lambda_param;
    double C_param, alpha_param;
    std::array<double, 9> B_param;
    std::array<Point<dim>, dim> aniso_fun_points;

    struct DEntry {
        int value;
        // expand with different function on D boundary
    };
    std::vector<DEntry> D_entries;
};

std::vector<Work> parse_file(const std::string& path);
}
