#include <deal.II/base/convergence_table.h>

#include "headers/NeoHooke.hpp"
#include "headers/Guccione.hpp"
#include "headers/TestConditions.hpp"
#include "headers/Work.hpp"

#include <filesystem>
#include <chrono>
using namespace std::chrono_literals;
namespace stdc = std::chrono;

void print_time(const stdc::nanoseconds& time, const ConditionalOStream& out);
void RunConvergenceTest();

int main(int argc, char *argv[])
{
    using namespace pde;
    constexpr unsigned int r = 1;

    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

    // This MPI process.
    const unsigned int mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));
    // Parallel output stream.
    const ConditionalOStream pcout(std::cout, mpi_rank == 0);

    pcout << argc << '\n';
    if (argc < 2)
    {
	    pcout << "Provide file with work\n";
	    pcout << "Usage: ./PDE-06 <path/to/config_file.txt>" << std::endl;
        return 1;
    }

    // if (argv[1][0] == '?')
    // {
    //     RunConvergenceTest();
    //     return 0;
    // }

    // -- PARSING WORKS --
    std::string filename(argv[1]);

    if (!std::filesystem::exists(filename))
    {
	    pcout << "File provided doesn't exist\n";
        return 1;
    }

    std::vector<Work> work = parse_file(filename, pcout);

    // -- OPENING TIME FILE --
    std::ofstream ftime;
    if (mpi_rank == 0)
        ftime.open("times.txt");
    
    const ConditionalOStream time_file(ftime, mpi_rank == 0);

    // -- RUNNING WORKS (JOBS) --
    pcout << "Work size " << work.size() << '\n';
    int failed = 0;
    for (auto& w : work)
    {
        pcout << "Starting work:";
        if (!w.name.empty()) pcout << " " << w.name << "\n";
        else pcout << "\n";
        if (w.passed) pcout << "\033[32m!THIS TEST WAS MARKED AS PASSED!\033[0m'\n";

	    // -- Initializing Mesh --
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
     
        pcout << "Solver iterations limit: " << w.iterations << '\n';

	    // -- DIRICHLET CONDITIONS --
        std::map<types::boundary_id, const Function<dim>*> dirichlet_conditions;
    	Functions::ZeroFunction<dim> zero_function(dim);
        TestDirichletConditions::SinXYFunction sin_function;

        pcout << "Diritto boundary on ids:";
        for (auto d : w.dirichlet_conditions)
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
        for (auto d : w.neumann_ids)
            pcout << ' ' << d;
        pcout << "\n";

	    // -- Setting up newton damping --
        pcout << "Newton scaling: " << w.newton_damping << '\n';

	    // -- NEUMANN CONDITIONS --
        std::function<Point<dim>(const Point<dim> &)> neumann_condition;

        try {
            if (w.neumann_func_param == "") {
                pcout << "No parameter for Neumann Condition found.\n" << "Falling back to standard: tau = 0.5" << std::endl;
                TestNeumannConditions::initialize();
            } else {
                TestNeumannConditions::initialize(std::stod(w.neumann_func_param));
            }
            neumann_condition = TestNeumannConditions::choose_neumann_function(w.neumann_func);
        } catch (std::invalid_argument &ia) {
            pcout << "Invalid Argument: " << ia.what() << std::endl;
        } catch (std::runtime_error &e) { 
            pcout << e.what() << " skipping work" << std::endl;
            continue;
        }
	
	    // -- FORCING TERM --
        pcout << "Forcing Term: " << w.forcing_term << '\n';

        MechanicalDisplacement::Config config{
                w.iterations, w.output_filename,
                std::move(mesh_src), r,
                w.newton_damping,
                neumann_condition, w.neumann_ids,
                dirichlet_conditions, 
		        TestForcingFunctions::choose_forcing_term(w.forcing_term)
            };

        pcout << "Time: " << w.time << '\n';
        int number_executions = w.time ? 10 : 1;

        try {
            if (w.material == Work::MaterialType::NeoHooke) // -- NEOHOOKE --
            {
                pcout << "NeoHooke Problem\n";
                
                Work::NeoHookeData params = std::get<Work::NeoHookeData>(w.problem_params);
                pcout << "mu = " << params.mu << ", lambda = " << params.lambda << '\n';

                NeoHooke problem = NeoHooke(std::move(config), pcout, mpi_rank, params.mu, params.lambda);
                
                stdc::nanoseconds total_time = 0ns;
                for (int ex = 0; ex < number_executions; ex++)
                {
                    auto start = stdc::high_resolution_clock::now();
                    problem.setup();
                    problem.solve();
                    auto end = stdc::high_resolution_clock::now();
                    total_time += end - start;
                    problem.output();
                }
                pcout << "Execution time: ";
                float avg = float(total_time.count()) / number_executions;
                stdc::nanoseconds avgns{(uint64_t)avg};
                print_time(avgns, pcout);
                pcout << '\n';

                time_file << "Time "; 
                print_time(avgns, time_file);
                time_file << '\n';
            }
            else // -- GUCCIONE --
            {
                pcout << "Guccione Problem\n";
                
                Work::GuccioneData params = std::get<Work::GuccioneData>(w.problem_params);
                pcout << "c = " << params.c << '\n';
                pcout << "b = ";
                for (unsigned int i = 0; i < params.b.size() - 1; i++)
                    pcout << params.b[i] << ", ";
                pcout << params.b[params.b.size() - 1] << '\n';

                const AnisotropicFunctionType aniso_fun = [&params](const Point<dim>&) {
                    return std::array<Point<dim>, dim>{ 
                        params.aniso_fun_points[0], params.aniso_fun_points[1], params.aniso_fun_points[2] 
                    };
                };

                Guccione problem = Guccione(std::move(config), pcout, mpi_rank, params.c, params.b, aniso_fun);
                
                stdc::nanoseconds total_time = 0ns;
                for (int ex = 0; ex < number_executions; ex++)
                {
                    auto start = stdc::high_resolution_clock::now();
                    problem.setup();
                    problem.solve();
                    auto end = stdc::high_resolution_clock::now();
                    total_time += end - start;
                    problem.output();
                }
                pcout << "Execution time: ";
                float avg = float(total_time.count()) / number_executions;
                stdc::nanoseconds avgns{(uint64_t)avg};
                print_time(avgns, pcout);
                pcout << '\n';

		        time_file << "Time "; 
                print_time(avgns, time_file);
                time_file << '\n';
            }
        } catch (const dealii::ExceptionBase& e) {
            failed++;
            pcout << "Deal.II exception raised:\n" << e.what() << "\nAborting!" << std::endl;
        } catch (std::exception& e) {
            failed++;
            pcout << "Exception raised:\n" << e.what() << "\nAborting!" << std::endl;
        }

        pcout << "\n\n\n\n";
    }

    if (failed != 0)
        pcout << "FAILED TASKS: " << failed << '\n';

    ftime.close();

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

    out << std::fixed << std::setprecision(3) << value << unit;

    out.get_stream().copyfmt(oldState);
}


namespace pde
{
class ExactSolution : public Function<dim>
{
public:
    ExactSolution() : Function<dim>(dim) {}

    virtual void vector_value(const Point<dim> &p, Vector<double> &values) const override
    {
        values[0] = p[2] * (p[2] - 1);
        values[1] = p[2] * (p[2] - 1);
        values[2] = 0;
    }

    virtual double value(const Point<dim> &p, const unsigned int component = 0) const override
    {
        if (component == 0) return p[2] * (p[2] - 1);
        if (component == 1) return p[2] * (p[2] - 1);
        return 0;
    }

    virtual void vector_gradient(const Point<dim> &p, std::vector<Tensor<1, dim>> &gradients) const override
    {
        gradients[0].clear();
        gradients[0][2] = 2 * p[2] - 1;
        gradients[1].clear();
        gradients[1][2] = 2 * p[2] - 1;
        gradients[2].clear();
    }
};
}

void RunConvergenceTest()
{
    using namespace pde;

    const unsigned int mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));
    const ConditionalOStream pcout(std::cout, mpi_rank == 0);

    ConvergenceTable table;

    const std::vector<unsigned int> refinement_values = {1, 2};

    const ExactSolution exact_solution;

    for (const auto &refinement : refinement_values)
    {
        std::unique_ptr<MeshGenerator<dim>> mesh_src;
        mesh_src = std::make_unique<CubeGenerator<dim>>(refinement);
        
        std::function<Point<dim>(const Point<dim> &)> neumann_condition = [] (const Point<dim> & p) { return p; }; // this is not used anyway
        std::unordered_set<int> N_ids = {};

        std::map<types::boundary_id, const Function<dim>*> dirichlet_conditions;
        Functions::ZeroFunction<dim> zero_function(dim);
        dirichlet_conditions[4] = &zero_function;
        dirichlet_conditions[5] = &zero_function;

        ForcingTermType forcing_term = [](const Point<dim> &) {
            Tensor<1, dim> res;
            res[0] = -2;
            res[1] = -2;
            res[2] = 0;	
            return res;
        };

        using namespace std::string_literals;

        MechanicalDisplacement::Config config{
                10000, "../../out/neo/cube/conv"s + std::to_string(refinement),
                std::move(mesh_src), 1, 0.3,
                neumann_condition, N_ids,
                dirichlet_conditions, 
                forcing_term
            };

        NeoHooke problem = NeoHooke(std::move(config), pcout, mpi_rank, 1, 2);
        
        problem.setup();
        problem.solve();
        problem.output();

        const double error_L2 = problem.compute_error(VectorTools::L2_norm, exact_solution);
        const double error_H1 = problem.compute_error(VectorTools::H1_norm, exact_solution);

        table.add_value("ref", refinement);
        table.add_value("L2", error_L2);
        table.add_value("H1", error_H1);
    }

    table.evaluate_all_convergence_rates(ConvergenceTable::reduction_rate_log2);

    table.set_scientific("L2", true);
    table.set_scientific("H1", true);

    if (mpi_rank == 0)
        table.write_text(std::cout);
}