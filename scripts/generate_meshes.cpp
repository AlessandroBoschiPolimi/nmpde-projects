#include "MeshGenerator.hpp"

#include <memory>

int main(int argc, char* argv[]) {

    if(argc < 5) {
	std::cout << "Usage: ./generate_meshes.cpp <[c]ube|[r]od> <[h]ex|[t]etra> <mesh_size_n> num_cells_msh_1 ... num_cells_msh_n <output_dir>" << std::endl;
	return 1;
    }

    using namespace pde;
    constexpr unsigned int dim = 3;

    char mesh_geom = *argv[1];
    char mesh_type = *argv[2];
    int mesh_size_num = std::stoi(argv[3]);
    std::vector<unsigned int> num_cells_per_mesh;
    for(int i = 4; i < mesh_size_num + 4 && i < argc; i++) {
	std::cout << argv[i] << std::endl;
	num_cells_per_mesh.push_back( static_cast<unsigned int>(std::stoi(argv[i])) );
    }
    
    std::string mesh_out_dir = "./mesh/";
    if(argc > mesh_size_num + 4) {
	mesh_out_dir = std::string(argv[mesh_size_num + 4]);
    }
    
    std::unique_ptr<MeshGenerator<dim>> generator;
    std::cout << "Setting Meshes (Type) " << mesh_type << std::endl;

    switch(mesh_type) {
	case 'h':
		std::cout << "Setting Meshes (Geom) " << mesh_geom << std::endl;
	    switch(mesh_geom) {
		case 'c':
		    generator = std::make_unique<CubeGenerator<dim>>(mesh_out_dir, num_cells_per_mesh);
		    break;
		case 'r':
		    generator = std::make_unique<RodGenerator<dim>>(mesh_out_dir, num_cells_per_mesh);
		    break;
		default:
		    std::cout << "ERROR: mesh geometry not supported. " << mesh_geom << std::endl;
	    }
	    break;
	case 't':
		std::cout << "Setting Meshes (Geom) " << mesh_geom << std::endl;
	    switch(mesh_geom) {
		// case 'e':
		//     generator = std::make_unique<MeshLoader<dim>>();
		//     break;
		case 'c':
		    generator = std::make_unique<CubeGeneratorWithSimplices<dim>>(mesh_out_dir, num_cells_per_mesh);
		    break;
		case 'r':
		    generator = std::make_unique<RodGeneratorWithSimplices<dim>>(mesh_out_dir, num_cells_per_mesh);	
		    break;
		default:
		    std::cout << "ERROR: mesh geometry not supported. " << mesh_geom << std::endl;
	    }
	    break;
	default:
	    std::cout << "ERROR: mesh type not supported " << mesh_type << std::endl;
    }
    
	std::cout << "Generating Mesh... " << std::endl;
    generator->Generate();
    return 0;
}
