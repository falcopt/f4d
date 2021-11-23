#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>

#define VERBOSE_LEVEL 1  // SPH verbose level

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
    auto rvnd2 = cobra::RandomizedVariableNeighborhoodDescent(instance, move_generators, {cobra::TLCH}, rand_engine, tolerance);

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

    /**
     * DIMACS config
     */

    int time_base;
    int time_fact;
    int runs;
    bool restart = false;
    int sph_timelimit;

    if (instance.get_customers_num() <= 1000) {
        time_base = 120;
        time_fact = 2;
        runs = 3;
        restart = true;
        sph_timelimit = 60;
    } else if (instance.get_customers_num() <= 10000) {
        time_base = 600;
        time_fact = 2;
        runs = 1;
        sph_timelimit = 360;
    } else {
        time_base = 1200;
        time_fact = 2;
        runs = 1;
        sph_timelimit = 720;
    }

    sph.set_timelimit(sph_timelimit);

    solution = best_solution;

    while (true) {

        for (int run = 0; run < runs; ++run) {

            if (restart) solution = cw_solution;

            solution = cos.coreopt(solution, time_base);
        }



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

        time_base *= time_fact;

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