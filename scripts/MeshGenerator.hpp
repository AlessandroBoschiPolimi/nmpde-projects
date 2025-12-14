#ifndef MESH_GENERATOR_HPP
#define MESH_GENERATOR_HPP

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/grid/grid_generator.h>

#include <filesystem>

namespace pde {

using namespace dealii;

enum Type
{
    Tetrahedra,
    Hexahedra
};

template <int dim>
class MeshGenerator
{
protected:
    const std::filesystem::path mesh_path;
    const std::vector<unsigned int> num_cells_per_mesh;
    void write_mesh(const Triangulation<dim>& mesh, const std::string& mesh_type, const unsigned int num_cells) const {
	const std::string mesh_file_name = mesh_path.parent_path().string() + "/" + "mesh-" + mesh_type + "-" + std::to_string(num_cells) + ".msh";
	GridOut           grid_out;
	std::ofstream     grid_out_file(mesh_file_name);
	grid_out.write_msh(mesh, grid_out_file);
	std::cout << "Mesh saved to " << mesh_file_name << std::endl;
	grid_out_file.close();
    };

public:
    MeshGenerator(const std::string &mesh_output_file_, const std::vector<unsigned int>& num_cells_per_mesh_): 
	mesh_path(mesh_output_file_), num_cells_per_mesh(num_cells_per_mesh_) {};
    virtual ~MeshGenerator() = default;

    virtual void Generate() const = 0;
    virtual Type ElementType() const = 0;
};

template <int dim>
class MeshLoader : public MeshGenerator<dim>
{
public:
    ~MeshLoader() override = default;

    void Generate() const override
    {
        Triangulation<dim> mesh;
        {
            GridIn<dim> grid_in;
            grid_in.attach_triangulation(mesh);

            std::ifstream grid_in_file("../mesh/mesh.msh");
            grid_in.read_msh(grid_in_file);
        }

        // // Then, we copy the triangulation into the parallel one.
        // {
        //     GridTools::partition_triangulation(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD), mesh);
        //     const auto construction_data = TriangulationDescription::Utilities::create_description_from_triangulation(mesh, MPI_COMM_WORLD);
        //     mesh.create_triangulation(construction_data);
        // }
    }

    Type ElementType() const override { return Type::Tetrahedra; }

private:
    static constexpr int outside_id = 1, top_id = 2, inside_id = 3;
};

template <int dim>
class RodGenerator : public MeshGenerator<dim>
{
public:
    using MeshGenerator<dim>::MeshGenerator;
    ~RodGenerator() override = default;

    void Generate() const override
    {
	for(const auto num_cells : this->num_cells_per_mesh) {
	    Triangulation<dim> mesh;
	    GridGenerator::subdivided_cylinder(mesh, num_cells, 0.1, 1.0);
	    this->write_mesh(mesh, "rod-hex", num_cells);
	}
    }

    Type ElementType() const override { return Type::Hexahedra; }

private:
    static constexpr int hull_id = 0, left_id = 1, right_id = 2;
};

template <int dim>
class RodGeneratorWithSimplices : public MeshGenerator<dim>
{
public:
    using MeshGenerator<dim>::MeshGenerator;
    ~RodGeneratorWithSimplices() override = default;

    void Generate() const override
    {
	for(const auto num_cells : this->num_cells_per_mesh) {
	    Triangulation<dim> mesh;
	    GridGenerator::subdivided_hyper_rectangle_with_simplices(
		mesh, {num_cells, num_cells, num_cells},
		Point<dim>(0,0,0), Point<dim>(2,1,1), true
	    );
	    this->write_mesh(mesh, "rod-tetra", num_cells);
	}
    }

    Type ElementType() const override { return Type::Tetrahedra; }
};

/* From the documentation:
    
Faces (quads in 3d): first the two faces with normals in x-, then y- and z-direction. 
For each two faces: first the face with normal in negative coordinate direction, then the one with normal in positive direction, 
i.e. the faces are ordered according to their normals pointing in -x, x, -y, y, -z, z direction.

Therefore, the faces are numbered in the ordering: left, right, front, back, bottom and top face:

*       *-------*        *-------*
*      /|       |       /       /|
*     / |   3   |      /   5   / |
*    /  |       |     /       /  |
*   *   |       |    *-------*   |
*   | 0 *-------*    |       | 1 *
*   |  /       /     |       |  /
*   | /   4   /      |   2   | /
*   |/       /       |       |/
*   *-------*        *-------*
* 
*/

template <int dim>
class CubeGenerator : public MeshGenerator<dim>
{
public:
    using MeshGenerator<dim>::MeshGenerator;
    ~CubeGenerator() override = default;

    void Generate() const override 
    {
	for(const auto num_cells : this->num_cells_per_mesh) {
	    Triangulation<dim> mesh;
	    GridGenerator::subdivided_hyper_cube(mesh, num_cells, 0.0, 1.0, true);
	    this->write_mesh(mesh, "cube-hex", num_cells);
	}
    }

    Type ElementType() const override { return Type::Hexahedra; }

};


template <int dim>
class CubeGeneratorWithSimplices : public MeshGenerator<dim>
{

public:
    using MeshGenerator<dim>::MeshGenerator; // Inheriting constructor
    ~CubeGeneratorWithSimplices() override = default;

    // Create the mesh using dealii generator, this gives us numbering of the faces
    void Generate() const override 
    {
	std::cout << "Generating Mesh for a Cube with Simplex support... " << std::endl;
	for(const auto num_cells : this->num_cells_per_mesh) {
		std::cout << "Initializing Mesh... " << num_cells << std::endl;
	    Triangulation<dim> mesh;
	    GridGenerator::subdivided_hyper_cube(mesh, num_cells, 0.0, 1.0, true);
	    this->write_mesh(mesh, "cube-tetra", num_cells);
	}
    }

    Type ElementType() const override { return Type::Tetrahedra; }
};

}

#endif // MESH_GENERATOR_HPP
