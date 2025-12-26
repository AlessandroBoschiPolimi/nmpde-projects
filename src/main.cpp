#include "headers/NeoHooke.hpp"
#include "headers/Guccione.hpp"
#include "headers/TestConditions.hpp"
#include "headers/Work.hpp"

#include <filesystem>
#include <chrono>
using namespace std::chrono_literals;
namespace stdc = std::chrono;

const pde::ForcingTermType select_forcing_term(std::string boolean_value) {
    using namespace pde::TestForcingFunctions;
    return std::stoi(boolean_value) ? bend_rod : null_forcing_term;
}

void print_time(const stdc::nanoseconds& time, const ConditionalOStream& out);

int main(int argc, char *argv[])
{
    using namespace pde;
    constexpr unsigned int r = 1;

    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

    // This MPI process.
    const unsigned int mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));
    // Parallel output stream.
    const ConditionalOStream pcout(std::cout, mpi_rank == 0);

    if (argc < 3)
    {
	    pcout << "Provide file with work\n";
	    pcout << "Usage: ./PDE-06 <path/to/config_file.txt> [FORCING_TERM <0|1>]" << std::endl;
        return 1;
    }

    std::string filename(argv[1]);

    if (!std::filesystem::exists(filename))
    {
	    pcout << "File provided doesn't exist\n";
        return 1;
    }

    std::vector<Work> work = parse_file(filename);

    pcout << "Work size " << work.size() << '\n';
    for (auto& w : work)
    {
        auto start = stdc::high_resolution_clock::now();

        pcout << "Starting work:\n";
        std::unique_ptr<MeshGenerator<dim>> mesh_src;
        if (w.geometry == Work::GeometryType::File) {
            pcout << "Using mesh from file: " << std::get<std::filesystem::path>(w.mesh_param.value()) << '\n';
            mesh_src = std::make_unique<MeshLoader<dim>>(std::get<std::filesystem::path>(w.mesh_param.value()));
        } else if (w.geometry == Work::GeometryType::Cube) {
            pcout << "Using cube mesh, refinement " << std::get<unsigned int>(w.mesh_param.value()) << '\n';
            mesh_src = std::make_unique<CubeGenerator<dim>>(std::get<unsigned int>(w.mesh_param.value()));
        } else if (w.geometry == Work::GeometryType::Rod) {
            pcout << "Using rod mesh\n";
            mesh_src = std::make_unique<RodGenerator<dim>>();
        }
     
        pcout << "Solver iterations limit " << w.iterations << '\n';

        std::map<types::boundary_id, const Function<dim>*> dirichlet_conditions;
    	Functions::ZeroFunction<dim> zero_function(dim);
        TestDirichletConditions::SinXYFunction<dim> sin_function;

        pcout << "Diritto boundary on ids:";
        for (auto d : w.D_entries)
        {
            pcout << ' ' << d.value << ' ' << d.function;
            if (d.function == "zero")
                dirichlet_conditions[d.value] = &zero_function;
            else if (d.function == "sin")
                dirichlet_conditions[d.value] = &sin_function;
            else {
                pcout << "Unknown Dirichlet boundary condition, skipping work\n";
                continue;
            }
        }
	    pcout << "\n";
        pcout << "UomoNuovo boundary on ids:";
        for (auto d : w.N_values)
            pcout << ' ' << d;
        pcout << "\n";
        if (w.newton_damping)
            pcout << "Newton scaling " << w.newton_scaling << '\n';


        std::function<Point<dim>(const Point<dim> &)> neumann_condition;

        try {
            if (w.N_data == "") {
                pcout << "No parameter for Neumann Condition found.\n"
                    << "Falling back to standard: tau = 0.5" << std::endl;
                TestNeumannConditions::initialize();
            } else {
                TestNeumannConditions::initialize(std::stod(w.N_data));
            }
            neumann_condition = TestNeumannConditions::choose_neumann_function(w.N_label);
        } catch(std::invalid_argument &ia) {
            pcout << "Invalid Argument: " << ia.what() << std::endl;
        } catch(std::runtime_error &e) { 
            pcout << e.what() << " skipping work" << std::endl;
            continue;
        }

        MechanicalDisplacement::Config config{
                w.iterations, w.output_filename, 
                std::move(mesh_src), r,
                w.newton_damping, w.newton_scaling,
                neumann_condition, w.N_values,
                dirichlet_conditions, select_forcing_term(argv[2])
            };

        // TODO: fix/add forcing term to work
        try {
            if (w.material == Work::MaterialType::NeoHooke)
            {
                pcout << "NeoHooke Problem\n";
                
                Work::NeoHookeData params = std::get<Work::NeoHookeData>(w.problem_params);
                pcout << "C = " << params.C << ", lambda = " << params.lambda << '\n';

                NeoHooke problem = NeoHooke(std::move(config), pcout, mpi_rank, params.C, params.lambda);
                problem.setup();
                problem.solve();
                problem.output();
            }
            else
            {
                pcout << "Guccione Problem\n";
                
                Work::GuccioneData params = std::get<Work::GuccioneData>(w.problem_params);
                pcout << "c = " << params.c << '\n';
                pcout << "b = ";
                for (unsigned int i = 0; i < params.b.size() - 1; i++)
                    pcout << params.b[i] << ", ";
                pcout << params.b[params.b.size() - 1] << '\n';

                const AnisotropicFunctionType aniso_fun = [&params](const Point<dim>&){ 
                        return std::array<Point<dim>, dim>{ params.aniso_fun_points[0], params.aniso_fun_points[1], params.aniso_fun_points[2] };
                    };

                Guccione problem = Guccione(std::move(config), pcout, mpi_rank, params.c, params.b, aniso_fun);
                problem.setup();
                problem.solve();
                problem.output();
            }
        } catch (const dealii::ExceptionBase& e) {
            pcout << "Deal.II exception raised:\n" << e.what() << "\nAborting!" << std::endl;
        } catch (std::exception& e) {
            pcout << "Exception raised:\n" << e.what() << "\nAborting!" << std::endl;
        }

        auto end = stdc::high_resolution_clock::now();
        auto dt = end - start;
        pcout << "Execution time: ";
        print_time(dt, pcout);

        pcout << "\n\n\n\n";

    }

    return 0;
}


void print_time(const stdc::nanoseconds& time, const ConditionalOStream& out)
{
    double value = static_cast<double>(time.count());
    const char* unit = "ns";

    if (value >= 1'000'000'000.0) {
        value /= 1'000'000'000.0;
        unit = "s";
    }
    else if (value >= 1'000'000.0) {
        value /= 1'000'000.0;
        unit = "ms";
    }
    else if (value >= 1'000.0) {
        value /= 1'000.0;
        unit = "us";
    }

    std::ios oldState(nullptr);
    oldState.copyfmt(out.get_stream());

    out << std::fixed << std::setprecision(2) << value << unit;

    out.get_stream().copyfmt(oldState);
}