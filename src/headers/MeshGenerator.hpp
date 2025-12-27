#pragma once

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/grid/grid_generator.h>

#include <fstream>
#include <filesystem>

using namespace dealii;

enum Type
{
    Tetrahedra,
    Hexahedra
};

template <int dim>
class MeshGenerator
{
public:
    virtual ~MeshGenerator() = default;

    virtual void Generate(parallel::fullydistributed::Triangulation<dim>& mesh) const = 0;
    virtual Type ElementType() const = 0;
};

template <int dim>
class MeshLoader : public MeshGenerator<dim>
{
public:
    MeshLoader(const std::filesystem::path& filename_)
        : filename(filename_) {}
    ~MeshLoader() override = default;

    void Generate(parallel::fullydistributed::Triangulation<dim>& mesh) const override
    {
        Triangulation<dim> mesh_serial;

        {
            GridIn<dim> grid_in;
            grid_in.attach_triangulation(mesh_serial);

            std::ifstream grid_in_file(filename);
            grid_in.read_msh(grid_in_file);
        }

        // Then, we copy the triangulation into the parallel one.
        {
            GridTools::partition_triangulation(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD), mesh_serial);
            const auto construction_data = TriangulationDescription::Utilities::create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
            mesh.create_triangulation(construction_data);
        }
    }

    Type ElementType() const override { return Type::Tetrahedra; }

public:
    static constexpr int outside_id = 1, top_id = 2, inside_id = 3;

    const std::filesystem::path filename;
};

template <int dim>
class CubeGenerator : public MeshGenerator<dim>
{
public:
    const unsigned int precision = false;

    CubeGenerator(unsigned int precision_ = 1) : precision(std::max(precision_, (unsigned int)1)) {}
    ~CubeGenerator() override = default;

    void Generate(parallel::fullydistributed::Triangulation<dim>& mesh) const override 
    {
        // Create the mesh using dealii generator, this gives us numbering of the faces
		/* From the documentation:
			
		Faces (quads in 3d): first the two faces with normals in x-, then y- and z-direction. For each two faces: first the face with normal in negative coordinate direction, then the one with normal in positive direction, i.e. the faces are ordered according to their normals pointing in -x, x, -y, y, -z, z direction.
		
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
		Triangulation<dim> mesh_serial;
		GridGenerator::subdivided_hyper_cube(mesh_serial, precision * 6, 0.0, 1.0, true);
        
		// Then, we copy the triangulation into the parallel one.
		{
			GridTools::partition_triangulation(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD), mesh_serial);
			const auto construction_data = TriangulationDescription::Utilities::create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
			mesh.create_triangulation(construction_data);
		}
    }

    Type ElementType() const override { return Type::Hexahedra; }

};

template <int dim>
class RodGenerator : public MeshGenerator<dim>
{
public:
    ~RodGenerator() override = default;

    void Generate(parallel::fullydistributed::Triangulation<dim>& mesh) const override
    {
        // Create the mesh using dealii generator
		Triangulation<dim> mesh_serial;
		GridGenerator::subdivided_hyper_rectangle(mesh_serial, std::vector<unsigned int>{20, 3, 3}, Point<dim>(0, 0, 0), Point<dim>(10, 1, 1), true);

		// Then, we copy the triangulation into the parallel one.
		{
			GridTools::partition_triangulation(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD), mesh_serial);
			const auto construction_data = TriangulationDescription::Utilities::create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
			mesh.create_triangulation(construction_data);
		}
    }

    Type ElementType() const override { return Type::Hexahedra; }

    //static constexpr int hull_id = 0, left_id = 1, right_id = 2;
};
