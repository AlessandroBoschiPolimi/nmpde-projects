#include <optional>
#include <vector>
#include <string>
#include <unordered_set>

namespace pde {

struct Work {
    enum MaterialType { NeoHooke, Guccione };
    MaterialType material;
    enum GeometryType { File, Cube, Rod };
    GeometryType geometry;
    std::optional<std::string> filename; // only for file
    std::string output_filename;

    std::unordered_set<int> N_values;
    std::string N_label;
    std::string N_data;

    struct DEntry {
        int value;
        // expand with different function on D boundary
    };
    std::vector<DEntry> D_entries;
};

std::vector<Work> parse_file(const std::string& path);
}
