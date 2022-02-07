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
    int route = 1;
    for (sph::idx_t j : sol) {
        fmt::print("Route #{}: ", route++);
        for (sph::idx_t i : inst.get_col(j)) fmt::print("{} ", i + 1);
        fmt::print("\n");
    }
    fmt::print("Cost {}\n", sol.get_cost());
    fflush(stdout);
}


auto main([[maybe_unused]] int argc, char* argv[]) -> int {

    auto param = Parameters(std::string(argv[1]));

    param.set(TOKEN_ROUND, argv[2]);

    int T = std::stoi(argv[3]);

    auto rand_engine = std::mt19937(param.get_seed());

    const auto parser_type = param.get_parser();

    const auto round_costs = param.needs_round_costs();

    auto maybe_instance = (round_costs ? cobra::Instance::make<true>(param.get_instance_path()) : cobra::Instance::make<false>(param.get_instance_path()));

    if (!maybe_instance) {
        std::cout << "Error while parsing the instance '" << param.get_instance_path() << "'.\n";
        exit(EXIT_FAILURE);
    }

    const auto instance = std::move(maybe_instance.value());

    auto k = param.get_sparsification_rule_neighbors();
    auto knn_view = cobra::KNeighborsMoveGeneratorsView(instance, k);

    auto views = std::vector<cobra::AbstractMoveGeneratorsView*>();
    views.push_back(&knn_view);

    auto move_generators = cobra::MoveGenerators(instance, views);

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


    const auto solution_history_size = param.get_solution_cache_size();
    auto solution = cobra::Solution(instance, std::min(instance.get_vertices_num(), solution_history_size));

    cobra::clarke_and_wright(instance, solution, param.get_cw_lambda(), param.get_cw_neighbors());
    solution.print_dimacs();

    if (instance.get_vertices_num() > 400) {  // Here some hardcoded magic numbers for the DIMACS challenge

        auto kmin = greedy_first_fit_decreasing(instance);
        if (kmin < solution.get_routes_num()) {

            const auto routemin_iterations = param.get_routemin_iterations();

            solution = routemin(instance, solution, rand_engine, move_generators, kmin, routemin_iterations, tolerance);
        }
    }

    sph::SPHeuristic sph(instance.get_vertices_num() - 1);
    sph.set_new_best_callback(print_sol);
    sph.set_max_routes(500'000U);
    sph.set_keepcol_strategy(sph::SPP);

    TimeBasedCoreOptSolver cos(instance, param, rand_engine, move_generators, local_search, sph);

    auto best_solution = solution;

    int PHASE1 = 0;
    int PHASE2 = 1;

    std::array<float, 2> phase_time;
    std::array<int, 2> rounds;
    std::array<float, 2> filo_time;
    std::array<int, 2> filo_runs;
    std::array<float, 2> sph_time;
    std::array<std::string, 2> sa_init_factor;
    std::array<std::string, 2> sa_final_factor;

    if (instance.get_vertices_num() <= 200) {  // small

        phase_time[PHASE1] = 1.0;

        rounds[PHASE1] = 3;
        filo_time[PHASE1] = 0.8;
        filo_runs[PHASE1] = 8;

        rounds[PHASE2] = 0;
        filo_time[PHASE2] = 0.0;
        filo_runs[PHASE2] = 0;

    } else if (instance.get_vertices_num() <= 400) {  // medium

        phase_time[PHASE1] = 1.0;

        rounds[PHASE1] = 6;
        filo_time[PHASE1] = 0.8;
        filo_runs[PHASE1] = 8;

        rounds[PHASE2] = 0;
        filo_time[PHASE2] = 0.0;
        filo_runs[PHASE2] = 0;

    } else if (instance.get_vertices_num() <= 1001) {  // large

        phase_time[PHASE1] = 0.25;

        rounds[PHASE1] = 3;
        filo_time[PHASE1] = 0.6;
        filo_runs[PHASE1] = 3;

        rounds[PHASE2] = 5;
        filo_time[PHASE2] = 0.78;
        filo_runs[PHASE2] = 1;

    } else if (instance.get_vertices_num() <= 11001) {  // x large

        phase_time[PHASE1] = 0.12;

        rounds[PHASE1] = 1;
        filo_time[PHASE1] = 0.72;
        filo_runs[PHASE1] = 1;

        rounds[PHASE2] = 4;
        filo_time[PHASE2] = 0.85;
        filo_runs[PHASE2] = 4;

    } else {  // xx large

        phase_time[PHASE1] = 0.17;

        rounds[PHASE1] = 1;
        filo_time[PHASE1] = 0.80;
        filo_runs[PHASE1] = 1;

        rounds[PHASE2] = 2;
        filo_time[PHASE2] = 0.92;
        filo_runs[PHASE2] = 2;
    }

    phase_time[PHASE2] = 1 - phase_time[PHASE1];
    sph_time[PHASE1] = 1 - filo_time[PHASE1];
    sph_time[PHASE2] = 1 - filo_time[PHASE2];
    sa_init_factor[PHASE1] = "0.1";
    sa_init_factor[PHASE2] = "0.01";
    sa_final_factor[PHASE1] = "0.001";
    sa_final_factor[PHASE2] = "0.0001";

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

            if (refined_solution.get_cost() < solution.get_cost()) {
                solution = refined_solution;
            }

            if (solution.get_cost() < best_solution.get_cost()) {
                solution.print_dimacs();
                best_solution = solution;
            }
        }
    }


    return EXIT_SUCCESS;
}