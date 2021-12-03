#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>

#define VERBOSE_LEVEL 1  // SPH verbose level

#include "Config.hpp"
#include "LinKernighan.hpp"
#include "RuinAndRecreate.hpp"
#include "TimeBasedCoreOptSolver.hpp"
#include "arg_parser.hpp"
#include "cobra/Instance.hpp"
#include "cobra/LocalSearch.hpp"
#include "cobra/PrettyPrinter.hpp"
#include "cobra/RouteSet.hpp"
#include "cobra/Solution.hpp"
#include "cobra/SolutionAlgorithms.hpp"
#include "greedy_bpp.hpp"
#include "routemin.hpp"

inline std::string get_basename(const std::string& path) {
    return {std::find_if(path.rbegin(), path.rend(),
                         [](char c) {
                             return c == '/';
                         })
                .base(),
            path.end()};
}

void print_sol(sph::Instance& inst, sph::GlobalSolution& sol) {
#ifdef DIMACS
    int route = 1;
    for (sph::idx_t j : sol) {
        fmt::print("Route #{}: ", route++);
        for (sph::idx_t i : inst.get_col(j)) fmt::print("{} ", i + 1);
        fmt::print("\n");
    }
    fmt::print("Cost {}\n", sol.get_cost());
    fflush(stdout);
#endif
}

void store_to_file(const cobra::Instance& instance, const cobra::Solution& solution, const std::string& path) {

    auto out_stream = std::ofstream(path);

    for (auto route = solution.get_first_route(), idx = 1; route != cobra::Solution::dummy_route; route = solution.get_next_route(route), idx++) {
        out_stream << "Route #" << idx << ":";
        for (auto customer = solution.get_first_customer(route); customer != instance.get_depot(); customer = solution.get_next_vertex(customer)) {
            out_stream << " " << std::to_string(customer);
        }
        out_stream << "\n";
    }

    out_stream << "Cost " << std::to_string(solution.get_cost());
}


auto main(int argc, char* argv[]) -> int {

#ifndef NDEBUG
    std::cout << "******************************\n";
    std::cout << "Probably running in DEBUG mode\n";
    std::cout << "******************************\n\n";
#endif


#ifdef VERBOSE
    std::cout << "FILO (Fast ILS Localized Optimization)\n";
    std::cout << "--- 2.0-dimacs\n";
    std::cout << "\n";
    std::cout << "Build options\n";
    std::cout << "--- Verbose output ENABLED\n";
    #ifdef GUI
    std::cout << "--- Graphical interface ENABLED\n";
    #else
    std::cout << "--- Graphical interface DISABLED\n";
    #endif
    std::cout << "\n";
#endif

    auto param = parse_command_line_arguments(argc, argv);

    int T = std::stoi(argv[3]);
    Config cfg(argv[4]);

    // int T = 7200;
    // Config cfg("/home/acco/git/f4d/configs/baseline.config");

    param.set(TOKEN_OUTPATH, cfg.get_string("filo_outpath_prefix") + argv[4] + "/");

    std::filesystem::create_directories(param.get_outpath());

    const auto outfile = param.get_outpath() + get_basename(param.get_instance_path()) + "_seed-" + std::to_string(param.get_seed());

    const auto global_time_begin = std::chrono::high_resolution_clock::now();

    auto rand_engine = std::mt19937(param.get_seed());

#ifdef VERBOSE
    auto partial_time_begin = std::chrono::high_resolution_clock::now();
    auto partial_time_end = std::chrono::high_resolution_clock::now();

    std::cout << "Pre-processing the instance.\n";
    partial_time_begin = std::chrono::high_resolution_clock::now();
#endif

    const auto parser_type = param.get_parser();

    const auto round_costs = param.needs_round_costs();

    auto maybe_instance = (round_costs ? cobra::Instance::make<true>(param.get_instance_path()) : cobra::Instance::make<false>(param.get_instance_path()));


    if (!maybe_instance) {
        std::cout << "Error while parsing the instance '" << param.get_instance_path() << "'.\n";
        exit(EXIT_FAILURE);
    }

    const auto instance = std::move(maybe_instance.value());

#ifdef VERBOSE
    partial_time_end = std::chrono::high_resolution_clock::now();
    std::cout << "Done in " << std::chrono::duration_cast<std::chrono::milliseconds>(partial_time_end - partial_time_begin).count() << " milliseconds.\n\n";
#endif


#ifdef VERBOSE
    std::cout << "Computing a greedy upper bound on the n. of routes.\n";
    partial_time_begin = std::chrono::high_resolution_clock::now();
#endif
    auto kmin = greedy_first_fit_decreasing(instance);
#ifdef VERBOSE
    partial_time_end = std::chrono::high_resolution_clock::now();
    std::cout << "Done in " << std::chrono::duration_cast<std::chrono::milliseconds>(partial_time_end - partial_time_begin).count() << " milliseconds.\n";
    std::cout << "Around " << kmin << " routes should do the job.\n\n";

    std::cout << "Setting up MOVEGENERATORS data structures.\n";
    partial_time_begin = std::chrono::high_resolution_clock::now();
#endif

    auto k = param.get_sparsification_rule_neighbors();
    auto knn_view = cobra::KNeighborsMoveGeneratorsView(instance, k);

    auto views = std::vector<cobra::AbstractMoveGeneratorsView*>();
    views.push_back(&knn_view);

    auto move_generators = cobra::MoveGenerators(instance, views);

#ifdef VERBOSE
    partial_time_end = std::chrono::high_resolution_clock::now();
    std::cout << "Done in " << std::chrono::duration_cast<std::chrono::milliseconds>(partial_time_end - partial_time_begin).count() << " milliseconds.\n";
    const auto tot_arcs = instance.get_vertices_num() * instance.get_vertices_num();
    const auto move_gen_num = move_generators.get_raw_vector().size();
    const auto move_gen_perc = 100.0f * static_cast<float>(move_gen_num) / static_cast<float>(tot_arcs);
    std::cout << "Using at most " << move_generators.get_raw_vector().size() << " move-generators out of " << tot_arcs << " total arcs ";
    std::cout << std::fixed;
    std::cout << std::setprecision(2);
    std::cout << "(approx. " << move_gen_perc << "%):\n";
    std::cout << std::setprecision(10);
    std::cout << std::defaultfloat;
    std::cout << std::setw(10);
    std::cout << knn_view.get_number_of_moves() << " k=" << k << " nearest-neighbors arcs\n";
    std::cout << "\n";

    std::cout << "Setting up LOCALSEARCH data structures.\n";
    partial_time_begin = std::chrono::high_resolution_clock::now();
#endif

    const auto tolerance = param.get_tolerance();
    auto rvnd0 = cobra::RandomizedVariableNeighborhoodDescent(
        instance, move_generators,
        {cobra::E11,   cobra::E10,  cobra::TAILS, cobra::SPLIT, cobra::RE22B, cobra::E22,  cobra::RE20,  cobra::RE21,  cobra::RE22S, cobra::E21, cobra::E20,
         cobra::TWOPT, cobra::RE30, cobra::E30,   cobra::RE33B, cobra::E33,   cobra::RE31, cobra::RE32B, cobra::RE33S, cobra::E31,   cobra::E32, cobra::RE32S},
        rand_engine, tolerance);
    auto rvnd1 = cobra::RandomizedVariableNeighborhoodDescent(instance, move_generators, {cobra::EJCH, cobra::TLCH, cobra::STCH}, rand_engine, tolerance);

    auto local_search = cobra::VariableNeighborhoodDescentComposer(tolerance);
    local_search.append(&rvnd0);
    local_search.append(&rvnd1);

#ifdef VERBOSE
    partial_time_end = std::chrono::high_resolution_clock::now();
    std::cout << "Done in " << std::chrono::duration_cast<std::chrono::milliseconds>(partial_time_end - partial_time_begin).count() << " milliseconds.\n\n";
#endif

    const auto solution_history_size = param.get_solution_cache_size();
    auto solution = cobra::Solution(instance, std::min(instance.get_vertices_num(), solution_history_size));


#ifdef VERBOSE
    std::cout << "Running CLARKE&WRIGHT to generate an initial solution.\n";
    partial_time_begin = std::chrono::high_resolution_clock::now();
#endif
    cobra::clarke_and_wright(instance, solution, param.get_cw_lambda(), param.get_cw_neighbors());
    auto cw_solution = solution;
#ifdef VERBOSE
    partial_time_end = std::chrono::high_resolution_clock::now();
    std::cout << "Done in " << std::chrono::duration_cast<std::chrono::milliseconds>(partial_time_end - partial_time_begin).count() << " milliseconds.\n";
#endif


#ifdef VERBOSE
    std::cout << "Initial solution: obj = " << solution.get_cost() << ", n. of routes = " << solution.get_routes_num() << ".\n\n";
#endif


    assert(solution.is_feasible());

    if (kmin < solution.get_routes_num()) {
        const auto routemin_iterations = param.get_routemin_iterations();
#ifdef VERBOSE
        std::cout << "Running ROUTEMIN heuristic for at most " << routemin_iterations << " iterations.\n";
        std::cout << "Starting solution: obj = " << solution.get_cost() << ", n. of routes = " << solution.get_routes_num() << ".\n";
        partial_time_begin = std::chrono::high_resolution_clock::now();
#endif
        solution = routemin(instance, solution, rand_engine, move_generators, kmin, routemin_iterations, tolerance);
#ifdef VERBOSE
        std::cout << "Final solution: obj = " << solution.get_cost() << ", n. routes = " << solution.get_routes_num() << "\n";
        partial_time_end = std::chrono::high_resolution_clock::now();
        std::cout << "Done in " << std::chrono::duration_cast<std::chrono::seconds>(partial_time_end - partial_time_begin).count() << " seconds.\n\n";
#endif
        assert(solution.is_feasible());
    }

    sph::SPHeuristic sph(instance.get_vertices_num() - 1);
    sph.set_new_best_callback(print_sol);
    sph.set_max_routes(500'000U);
    sph.set_keepcol_strategy(sph::SPP);

    TimeBasedCoreOptSolver cos(instance, param, rand_engine, move_generators, local_search, sph);

    const auto fastopt_iterations = param.get_fastopt_iterations();
    auto best_solution = cos.fastopt(solution, fastopt_iterations);
#ifdef VERBOSE
    std::cout << "\n";
    std::cout << "Final solution found:\n";
    std::cout << "obj = " << solution.get_cost() << ", n. routes = " << solution.get_routes_num() << "\n\n";
#endif

    int small_size = cfg.get_int("small_size");
    int medium_size = cfg.get_int("medium_size");
    int large_size = cfg.get_int("large_size");
    int xlarge_size = cfg.get_int("xlarge_size");

    int PHASE1 = 0;
    int PHASE2 = 1;

    std::array<float, 2> phase_time;
    std::array<int, 2> rounds;
    std::array<float, 2> filo_time;
    std::array<int, 2> filo_runs;
    std::array<float, 2> sph_time;
    std::array<std::string, 2> sa_init_factor;
    std::array<std::string, 2> sa_final_factor;

    std::string size_prefix;
    if (instance.get_vertices_num() <= small_size) {
        size_prefix = "small";
    } else if (instance.get_vertices_num() <= medium_size) {
        size_prefix = "medium";
    } else if (instance.get_vertices_num() <= large_size) {
        size_prefix = "large";
    } else if (instance.get_vertices_num() <= xlarge_size) {
        size_prefix = "xlarge";
    } else {
        size_prefix = "xxlarge";
    }

    phase_time[PHASE1] = cfg.get_real(size_prefix + "_phase1_fract");
    phase_time[PHASE2] = 1 - phase_time[PHASE1];
    rounds[PHASE1] = cfg.get_int(size_prefix + "_phase1_rounds");
    rounds[PHASE2] = cfg.get_int(size_prefix + "_phase2_rounds");
    filo_time[PHASE1] = cfg.get_real(size_prefix + "_phase1_filo_time");
    filo_time[PHASE2] = cfg.get_real(size_prefix + "_phase2_filo_time");
    filo_runs[PHASE1] = cfg.get_int(size_prefix + "_phase1_filo_runs");
    filo_runs[PHASE2] = cfg.get_int(size_prefix + "_phase2_filo_runs");
    sph_time[PHASE1] = 1 - filo_time[PHASE1];
    sph_time[PHASE2] = 1 - filo_time[PHASE2];
    sa_init_factor[PHASE1] = cfg.get_string(size_prefix + "_phase1_sa_init_factor");
    sa_init_factor[PHASE2] = cfg.get_string(size_prefix + "_phase2_sa_init_factor");
    sa_final_factor[PHASE1] = cfg.get_string(size_prefix + "_phase1_sa_final_factor");
    sa_final_factor[PHASE2] = cfg.get_string(size_prefix + "_phase2_sa_final_factor");
    solution = best_solution;

    for (int phase = 0; phase < 2; ++phase) {

        int t_phase = T * phase_time[phase];

        if (t_phase == 0) continue;

        int t_filo = t_phase * filo_time[phase];
        int t_sph = (t_phase - t_filo) / rounds[phase];

        param.set(TOKEN_SA_INIT_FACTOR, sa_init_factor[phase]);
        param.set(TOKEN_SA_FINAL_FACTOR, sa_final_factor[phase]);

        int t_coreopt = t_filo / (filo_runs[phase] * rounds[phase]);

        for (int round = 0; round < rounds[phase]; ++round) {

            for (int run = 0; run < filo_runs[phase]; ++run) {

                solution = cos.coreopt(solution, t_coreopt);
                if (solution.get_cost() < best_solution.get_cost()) {
                    best_solution = solution;
                }
            }

            sph.set_timelimit(t_sph);
            std::vector<sph::idx_t> columns = sph.solve();

            cobra::Solution refined_solution(instance);
            refined_solution.reset();
            for (sph::idx_t col_idx : columns) {
                const sph::Column& col = sph.get_col(col_idx);
                int route = refined_solution.build_one_customer_route(col.front() + 1);
                for (size_t n = 1; n < col.size(); ++n) {
                    int customer = col[n] + 1;
                    refined_solution.insert_vertex_before(route, instance.get_depot(), customer);
                }
            }
            refined_solution.print_dimacs();

            if (refined_solution.get_cost() < solution.get_cost()) {
                solution = refined_solution;
            }

            if (solution.get_cost() < best_solution.get_cost()) {
                best_solution = solution;
                auto out_stream = std::ofstream(outfile + ".out");
                out_stream << std::setprecision(10);
                out_stream << best_solution.get_cost() << "\n";
                store_to_file(instance, best_solution, outfile + ".sol");
            }
        }
    }


    const auto global_time_end = std::chrono::high_resolution_clock::now();

#ifdef VERBOSE
    std::cout << "\n";
    std::cout << "Best solution found:\n";
    std::cout << "obj = " << best_solution.get_cost() << ", n. routes = " << best_solution.get_routes_num() << "\n";

    std::cout << "\n";
    std::cout << "Run completed in " << std::chrono::duration_cast<std::chrono::seconds>(global_time_end - global_time_begin).count() << " seconds ";
    std::cout << "(" << std::chrono::duration_cast<std::chrono::milliseconds>(global_time_end - global_time_begin).count() << " milliseconds).\n";
#endif

    auto out_stream = std::ofstream(outfile + ".out");
    out_stream << std::setprecision(10);
    out_stream << best_solution.get_cost() << "\t" << std::chrono::duration_cast<std::chrono::seconds>(global_time_end - global_time_begin).count() << "\n";
    store_to_file(instance, best_solution, outfile + ".sol");


    return EXIT_SUCCESS;
}