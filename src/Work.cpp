#include "headers/Work.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

namespace pde {

struct ConfigOptions {
    ConfigOptions():skip(false),time(false),passed(false),name("") {};
    bool skip;
    bool time;
    bool passed;
    std::string name;
};

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

std::vector<std::string> split(std::string s, const std::string& delimiter) {
	std::vector<std::string> tokens;
	size_t pos = 0;

	std::string token;
	while ((pos = s.find(delimiter)) != std::string::npos) {
		token = s.substr(0, pos);
		tokens.push_back(token);
		s.erase(0, pos + delimiter.length());
	}

	tokens.push_back(s);
	return tokens;
}

void parse_options(ConfigOptions& opt,
		   std::string str)
{
	std::vector<std::string> toks = split(str, " ");
	if (toks.empty()) toks.push_back(str); // Maybe there's only one option
	for (std::string &tok : toks) {
		if (tok == "@skip") {
			opt.skip = true;
			continue;
		}
		if (tok == "@time") {
			opt.time = true;
			continue;
		}
		if (tok.substr(0,5) == "@name") {
		    opt.name = tok.substr(tok.find("(") + 1);
		    opt.name.erase(opt.name.find(")"), std::string::npos);
		    continue;
		}
		if (tok == "@passed") {
		    opt.passed = true;
		    continue;
		}
	}
}

void apply_options(Work &current, ConfigOptions &opt) {
    current.time = opt.time;
    current.name = opt.name;
    current.passed = opt.passed;
}

std::vector<Work> parse_file(const std::filesystem::path& path
			     , const ConditionalOStream &pcout) {
	std::ifstream in(path);

	if (!in)
		throw std::runtime_error("Cannot open file: " + path.string());

	std::vector<Work> sections;
	std::string line;

	auto next_line = [&]() -> std::string {
		if (!std::getline(in, line))
			throw std::runtime_error("Unexpected end of file");
		return trim(line);
	};
	uint_fast16_t work_count = 1;

	while (std::getline(in, line)) {
		line = trim(line);
		if (line.empty())
			continue;

		if (line[0] == '#')
			continue;

		ConfigOptions opt;
		while (line[0] == '@') {
			parse_options(opt, line);
			line = next_line();
		}

		if (line.substr(0,5) != "-----")
			throw std::runtime_error("Expected section separator '-----'");

		Work sec;
		apply_options(sec, opt);
		bool skip = opt.skip;

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
				std::string filename;
				ss >> filename;
				if (filename.empty())
					throw std::runtime_error("Missing filename after 'file'");
				sec.mesh_param = filename;
			}
			else if (type == "cube") {
				sec.geometry = Work::GeometryType::Cube;
				unsigned int refinement = 1;
				ss >> refinement; // if there is no int at this point in the string: refinement is not modified
				sec.mesh_param = refinement;
			} 
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

		/* it line */
		line = next_line();
		{
			std::stringstream ss(line);
			std::string word;
			ss >> word;

			if (word == "it")
				ss >> sec.iterations;
			else throw std::runtime_error("Provide iteration count as fourth line");
		}

		/* scaling line */
		bool no_read = false;
		line = next_line();
		{
			std::stringstream ss(line);
			std::string word;
			ss >> word;

			if (word == "new_damn") {
				sec.newton_damping = true;
				ss >> sec.newton_scaling;
			}
			else
				no_read = true;
		}

		/* f line */
		if (!no_read)
			line = next_line();
		no_read = false;
		{
			std::stringstream ss(line);
			std::string word;
			ss >> word;

			if (word == "f") {
				ss >> sec.forcing_term;
			}
			else
				no_read = true;
		}

		/* N line */
		if (!no_read)
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
			if (!std::getline(in, line))
				break;

			line = trim(line);
			if (line.empty())
				continue;

			if (sec.material == Work::MaterialType::Guccione && line[0] == 'c')
			   break;

			if (sec.material == Work::MaterialType::NeoHooke && line[0] == 'C')
				break;

			if (line[0] != 'D')
				throw std::runtime_error("Expected D line");

			std::stringstream ss(line.substr(1));
			int value;
			std::string function;

			ss >> value >> function;
			if (function.empty())
				function = "zero";

			sec.D_entries.push_back({value, function});
		}
	
		// if in here already read line
		
		std::vector<std::string> toks;
		/* NeoHooke additional parameters */
		if (sec.material == Work::MaterialType::NeoHooke) {
			Work::NeoHookeData data;

			/* c param */
			{
				if (line.empty() || line[0] != 'C')
					throw std::runtime_error("Expected C line");
				std::string toks = line.substr(2);
				try {
					data.C = std::stod(toks);
				} catch (std::invalid_argument &invarg) {
					std::string err_msg = "[ERROR WHILE PARSING] at work " + std::to_string(work_count) + " while parsing c line " + invarg.what();
					throw std::invalid_argument(err_msg);
				}
			}

			line = next_line();
			/* lambda param */
			{
				if (line.empty() || line[0] != 'l')
					throw std::runtime_error("Expected lambda line");
				std::string toks = line.substr(2);
				try {
					data.lambda = std::stod(toks);
				} catch (std::invalid_argument &invarg) {
					std::string err_msg = "[ERROR WHILE PARSING] at work " + std::to_string(work_count) + " while parsing c line " + invarg.what();
					throw std::invalid_argument(err_msg);
				}
			}

			sec.problem_params = data;
		}
	
		/* Guccione additional parameters*/
		if(sec.material == Work::MaterialType::Guccione) {
			Work::GuccioneData data;

			/* c param */
			{
				if (line.empty() || line[0] != 'c')
					throw std::runtime_error("Expected c line");
				std::string toks = line.substr(2);
				try {
					data.c = std::stod(toks);
				} catch (std::invalid_argument &invarg) {
					std::string err_msg = "[ERROR WHILE PARSING] at work " + std::to_string(work_count) + " while parsing c line '" + toks + "': " + invarg.what();
					throw std::invalid_argument(err_msg);
				}
			}

			line = next_line();
			/* b param */
			{
				if (line.empty() || line[0] != 'b')
					throw std::runtime_error("Expected b line at work " + std::to_string(work_count));

				toks = split(line.substr(2), " ");
				for (size_t i = 0; i < toks.size(); i++) {
					try {
						data.b[i] = std::stod(toks[i]);
					} catch (std::invalid_argument &invarg) {
						std::string err_msg = "[ERROR WHILE PARSING] at work " + std::to_string(work_count) + " while parsing b line " + invarg.what();
						throw std::invalid_argument(err_msg);
					}
				}
			}

			line = next_line();
			/* aniso fun */
			{
				if (line.empty())
					throw std::runtime_error("Expected anisotropic function line at work " + std::to_string(work_count));

				toks = split(line, " ");

				if(toks[0] != "anfun") 
					throw std::runtime_error("Expected anisotropic function line at work " + std::to_string(work_count));

				std::vector<std::string> point_val;
				for (int i = 1; i < 4; i++) {
					point_val = split(toks[i], ",");
					data.aniso_fun_points[i-1] = Point<dim>(
						std::stod(point_val[0]), 
						std::stod(point_val[1]), 
						std::stod(point_val[2])
					);
				}
			}

			sec.problem_params = data;
		}

		if (!skip)
		    sections.push_back(std::move(sec));
		else
		    pcout << "Skipping " << opt.name << std::endl;

		work_count++;
	}

	return sections;
}

}
