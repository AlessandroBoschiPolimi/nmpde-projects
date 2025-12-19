#include "headers/NeoHooke.hpp"
#include "headers/Guccione.hpp"
#include "headers/TestConditions.hpp"
#include "headers/Work.hpp"

#include <filesystem>

const pde::ForcingTermType select_forcing_term(std::string boolean_value) {
    using namespace pde::TestForcingFunctions;
    return std::stoi(boolean_value) ? bend_rod : null_forcing_term;
}


int main(int argc, char *argv[])
{
    using namespace pde;
    constexpr unsigned int r = 1;


    if (argc < 3)
    {
	std::cout << "Provide file with work\n";
	std::cout << "Usage: ./PDE-06 <path/to/config_file.txt> [FORCING_TERM <0|1>]" << std::endl;
        return 1;
    }

    std::string filename(argv[1]);

    if (!std::filesystem::exists(filename))
    {
	std::cout << "File provided doesn't exist\n";
        return 1;
    }

    std::vector<Work> work = parse_file(filename);

    std::cout << "Work size " << work.size() << '\n';

    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

    // This MPI process.
    const unsigned int mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));
    // Parallel output stream.
    const ConditionalOStream pcout(std::cout, mpi_rank == 0);


    for (auto& w : work)
    {
        pcout << "Starting work:\n";
        std::unique_ptr<MeshGenerator<dim>> mesh_src;
        if (w.geometry == Work::GeometryType::File)
        {
            pcout << "Using mesh from file: " << w.filename.value() << '\n';
            mesh_src = std::make_unique<MeshLoader<dim>>(w.filename.value());
        }
        if (w.geometry == Work::GeometryType::Cube)
        {
            pcout << "Using cube mesh\n";
            mesh_src = std::make_unique<CubeGenerator<dim>>();
        }
        if (w.geometry == Work::GeometryType::Rod)
        {
            pcout << "Using rod mesh\n";
            mesh_src = std::make_unique<RodGenerator<dim>>();
        }
     
        pcout << "Solver iterations limit " << w.iterations << '\n';

        std::map<types::boundary_id, const Function<dim>*> boundary_functions;
    	Functions::ZeroFunction<dim> zero_function(dim);

        pcout << "Diritto boundary on ids:";
        for (auto d : w.D_entries)
        {
            pcout << ' ' << d.value;
            boundary_functions[d.value] = &zero_function;
        }
	    pcout << "\n";
        pcout << "UomoNuovo boundary on ids:";
        for (auto d : w.N_values)
            pcout << ' ' << d;
        pcout << "\n";


        std::function<Point<dim>(const Point<dim> &)> h;

        try {
            if(w.N_data == "") {
                pcout << "No parameter for Neumann Condition found.\n"
                    << "Falling back to standard: tau = 0.5" << std::endl;
                TestNeumannConditions::initialize();
            } else {
                TestNeumannConditions::initialize(std::stod(w.N_data));
            }
            h = TestNeumannConditions::choose_neumann_function(w.N_label);
        } catch(std::invalid_argument &ia) {
            pcout << "Invalid Argument: " << ia.what() << std::endl;
        } catch(std::runtime_error &e) { 
            pcout << e.what() << " skipping work" << std::endl;
            continue;
        } 

        // TODO: fix/add forcing term to work
        if (w.material == Work::MaterialType::NeoHooke)
        {
            pcout << "NeoHooke Problem\n";
            
            Work::NeoHookeData params = std::get<Work::NeoHookeData>(w.problem_params);
            pcout << "C = " << params.C << ", lambda = " << params.lambda << '\n';

            NeoHooke problem = NeoHooke(
                std::move(mesh_src), r, 
                boundary_functions, h, 
                w.N_values, select_forcing_term(argv[2]),
                w.output_filename, w.iterations,
                pcout, mpi_rank,
                params.C, params.lambda
            );
            problem.setup();
            problem.solve();
            problem.output();

            pcout << "\n\n\n\n";
        }
        else
        {
            pcout << "Guccione Problem\n";
            
            Work::GuccioneData params = std::get<Work::GuccioneData>(w.problem_params);
            pcout << "c = " << params.c << ", alpha " << params.alpha << '\n';
            pcout << "b = ";
            for (int i = 0; i < params.b.size() - 1; i++)
                pcout << params.b[i] << ", ";
            pcout << params.b[params.b.size() - 1] << '\n';

            const AnisotropicFunctionType aniso_fun = [&params](const Point<dim>&){ 
                    return std::array<Point<dim>, dim>(
                        { params.aniso_fun_points[0], params.aniso_fun_points[1], params.aniso_fun_points[2] }
                    );
                };

            Guccione problem = Guccione(
                std::move(mesh_src), r,
                boundary_functions, h,
                w.N_values, select_forcing_term(argv[2]),
                w.output_filename, w.iterations,
                pcout, mpi_rank,
                params.c, params.b,
                aniso_fun, params.alpha
            );
            problem.setup();
            problem.solve();
            problem.output();

            pcout << "\n\n\n\n";
        }
    }

    return 0;
}
