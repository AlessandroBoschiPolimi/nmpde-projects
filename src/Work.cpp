#include "headers/Work.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

namespace pde {

std::string trim(const std::string& s) {
    auto start = s.find_first_not_of(" \t\r\n");
    auto end   = s.find_last_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    return s.substr(start, end - start + 1);
}

std::unordered_set<int> parse_csv_ints(const std::string& s) {
    std::unordered_set<int> out;
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, ',')) {
        out.insert(std::stoi(trim(token)));
    }
    return out;
}

std::vector<Work> parse_file(const std::string& path) {
    std::ifstream in(path);

    if (!in)
        throw std::runtime_error("Cannot open file: " + path);

    std::vector<Work> sections;
    std::string line;

    auto next_line = [&]() -> std::string {
        if (!std::getline(in, line))
            throw std::runtime_error("Unexpected end of file");
        return trim(line);
    };

    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty())
            continue;

        if (line != "-----")
            throw std::runtime_error("Expected section separator '-----'");

        Work sec;

        /* Material */
        line = next_line();
        if (line == "NeoHooke") sec.material = Work::MaterialType::NeoHooke;
        else if (line == "Guccione") sec.material = Work::MaterialType::Guccione;
        else throw std::runtime_error("Unknown material: " + line);

        /* Geometry */
        line = next_line();
        {
            std::stringstream ss(line);
            std::string type;
            ss >> type;

            if (type == "file") {
                sec.geometry = Work::GeometryType::File;
                ss >> sec.filename.emplace();
                if (!sec.filename || sec.filename->empty())
                    throw std::runtime_error("Missing filename after 'file'");
            }
            else if (type == "cube") sec.geometry = Work::GeometryType::Cube;
            else if (type == "rod")  sec.geometry = Work::GeometryType::Rod;
            else throw std::runtime_error("Unknown geometry: " + type);
        }

        /* Output file */
        line = next_line();
        {
            std::stringstream ss(line);
            std::string word;
            ss >> word;

            if (word == "out") {
                ss >> sec.output_filename;
                if (sec.output_filename.empty())
                    throw std::runtime_error("Missing filename after 'out'");
            }
            else throw std::runtime_error("Provide output filename as third line");
        }

        /* N line */
        line = next_line();
        if (line.empty() || line[0] != 'N')
            throw std::runtime_error("Expected N line");

        sec.N_values = parse_csv_ints(trim(line.substr(1)));

        /* N data line */
        line = next_line();
        {
            std::stringstream ss(line);
            ss >> sec.N_label;
            std::getline(ss, sec.N_data);
            sec.N_data = trim(sec.N_data);
        }

        /* D lines */
        while (true) {
            std::streampos pos = in.tellg();
            if (!std::getline(in, line))
                break;

            line = trim(line);
            if (line.empty())
                continue;

            if (line == "-----") {
                in.seekg(pos);
                break;
            }

            if (line[0] != 'D')
                throw std::runtime_error("Expected D line");

            std::stringstream ss(line.substr(1));
            int value;
            std::string zero;

            ss >> value >> zero;
            if (zero != "zero")
                throw std::runtime_error("Expected 'zero' after D value");

            sec.D_entries.push_back({value});
        }

        sections.push_back(std::move(sec));
    }

    return sections;
}

}
