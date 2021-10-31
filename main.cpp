#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>

#include "LinKernighan.hpp"
#include "RuinAndRecreate.hpp"
#include "arg_parser.hpp"
#include "cobra/Instance.hpp"
#include "cobra/LocalSearch.hpp"
#include "cobra/NeighborAcceptance.hpp"
#include "cobra/PrettyPrinter.hpp"
#include "cobra/RouteSet.hpp"
#include "cobra/Solution.hpp"
#include "cobra/SolutionAlgorithms.hpp"
#include "cobra/Welford.hpp"
#include "greedy_bpp.hpp"
#include "routemin.hpp"

#ifdef GUI
    #include "Renderer.hpp"
#endif

inline std::string get_basename(const std::string& path) {
    return {std::find_if(path.rbegin(), path.rend(),
                         [](char c) {
                             return c == '/';
                         })
                .base(),
            path.end()};
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

void print_progress(cobra::PrettyPrinter& printer, const cobra::Instance& instance, int iter, const int iterations,
                    std::chrono::time_point<std::chrono::high_resolution_clock>& main_opt_loop_begin_time, std::vector<float>& gamma, std::vector<int>& omega,
                    cobra::SimulatedAnnealing& sa, cobra::Solution& best_solution, std::optional<float> bks) {


    const auto progress = 100.0f * static_cast<float>(iter + 1) / static_cast<float>(iterations);
    const auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - main_opt_loop_begin_time).count();
    const auto iter_per_second = static_cast<float>(iter + 1) / (static_cast<float>(elapsed_seconds) + 0.01f);
    const auto remaining_iter = iterations - iter;
    const auto estimated_rem_time = static_cast<float>(remaining_iter) / iter_per_second;

    auto gamma_mean = 0.0f;
    for (auto i = instance.get_vertices_begin(); i < instance.get_vertices_end(); i++) {
        gamma_mean += gamma[i];
    }
    gamma_mean = (gamma_mean / static_cast<float>(instance.get_vertices_num()));

    auto omega_mean = 0.0f;
    for (auto i = instance.get_customers_begin(); i < instance.get_customers_end(); i++) {
        omega_mean += omega[i];
    }
    omega_mean /= static_cast<float>(instance.get_customers_num());


    printer.print(progress, iter + 1, bks.has_value() ? 100.0f * (best_solution.get_cost() - *bks) / *bks : best_solution.get_cost(),
                  best_solution.get_routes_num(), iter_per_second, estimated_rem_time, gamma_mean, omega_mean, sa.get_temperature());
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

    auto maybe_instance = (round_costs ? cobra::Instance::make<cobra::TSPLIB95, true>(param.get_instance_path())
                                       : cobra::Instance::make<cobra::TSPLIB95, false>(param.get_instance_path()));


    if (!maybe_instance) {
        std::cout << "Error while parsing the instance '" << param.get_instance_path() << "'.\n";
        exit(EXIT_FAILURE);
    }

    const auto instance = std::move(maybe_instance.value());

#ifdef VERBOSE
    partial_time_end = std::chrono::high_resolution_clock::now();
    std::cout << "Done in " << std::chrono::duration_cast<std::chrono::milliseconds>(partial_time_end - partial_time_begin).count() << " milliseconds.\n\n";

    std::cout << "Computing mean arc cost.\n";
    partial_time_begin = std::chrono::high_resolution_clock::now();
#endif

    auto mean_arc_cost = 0.0;
    for (auto i = instance.get_vertices_begin(); i < instance.get_vertices_end() - 1; i++) {
        for (auto j = i + 1; j < instance.get_vertices_end(); j++) {
            mean_arc_cost += instance.get_cost(i, j);
        }
    }
    mean_arc_cost /= (instance.get_vertices_num() * (instance.get_vertices_num() - 1) / 2.0);

#ifdef VERBOSE
    partial_time_end = std::chrono::high_resolution_clock::now();
    std::cout << "Done in " << std::chrono::duration_cast<std::chrono::milliseconds>(partial_time_end - partial_time_begin).count() << " milliseconds.\n";
    std::cout << "Mean arc cost value is " << mean_arc_cost << ".\n";

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
        solution = routemin(instance, solution, rand_engine, move_generators, kmin, routemin_iterations, tolerance, round_costs);
#ifdef VERBOSE
        std::cout << "Final solution: obj = " << solution.get_cost() << ", n. routes = " << solution.get_routes_num() << "\n";
        partial_time_end = std::chrono::high_resolution_clock::now();
        std::cout << "Done in " << std::chrono::duration_cast<std::chrono::seconds>(partial_time_end - partial_time_begin).count() << " seconds.\n\n";
#endif
        assert(solution.is_feasible());
    }

    const auto coreopt_iterations = param.get_coreopt_iterations();

    auto best_solution = solution;

    const auto gamma_base = param.get_gamma_base();
    auto gamma = std::vector<float>(instance.get_vertices_num(), gamma_base);
    auto gamma_counter = std::vector<int>(instance.get_vertices_num(), 0);

    const auto delta = param.get_delta();
    auto average_number_of_vertices_accessed = cobra::Welford();

    auto gamma_vertices = std::vector<int>();
    for (auto i = instance.get_vertices_begin(); i < instance.get_vertices_end(); i++) {
        gamma_vertices.emplace_back(i);
    }
    move_generators.set_active_percentage(gamma, gamma_vertices);

    auto ruined_customers = std::vector<int>();

#ifdef VERBOSE
    std::cout << "Running COREOPT for " << coreopt_iterations << " iterations.\n";

    auto welford_rac_before_shaking = cobra::Welford();
    auto welford_rac_after_shaking = cobra::Welford();
    auto welford_local_optima = cobra::Welford();
    auto welford_shaken_solutions = cobra::Welford();
    auto printer = cobra::PrettyPrinter({
        {"%", cobra::PrettyPrinter::Field::Type::REAL, 5, " "},
        {"Iterations", cobra::PrettyPrinter::Field::Type::INTEGER, 10, " "},
        {param.get_bks_value().has_value() ? "Gap" : "Objective",
         round_costs && !param.get_bks_value().has_value() ? cobra::PrettyPrinter::Field::Type::INTEGER : cobra::PrettyPrinter::Field::Type::REAL, 10, " ", 4},
        {"Routes", cobra::PrettyPrinter::Field::Type::INTEGER, 6, " "},
        {"Iter/s", cobra::PrettyPrinter::Field::Type::REAL, 10, " "},
        {"Eta (s)", cobra::PrettyPrinter::Field::Type::REAL, 10, " "},
        {"Gamma", cobra::PrettyPrinter::Field::Type::REAL, 5, " "},
        {"Omega", cobra::PrettyPrinter::Field::Type::REAL, 6, " "},
        {"Temp", cobra::PrettyPrinter::Field::Type::REAL, 6, " "},
    #ifdef LO_STORAGE
        {"Solutions", cobra::PrettyPrinter::Field::Type::INTEGER, 10, " "},
        {"Routes", cobra::PrettyPrinter::Field::Type::INTEGER, 10, " "},
    #endif
    });

    auto main_opt_loop_begin_time = std::chrono::high_resolution_clock::now();

    auto elapsed_minutes = 0;
#endif

    auto rr = RuinAndRecreate(instance, rand_engine);

    const auto intensification_lb = param.get_shaking_lb_factor();
    const auto intensification_ub = param.get_shaking_ub_factor();

    const auto mean_solution_arc_cost = solution.get_cost() /
                                        (static_cast<float>(instance.get_customers_num()) + 2.0f * static_cast<float>(solution.get_routes_num()));

    auto shaking_lb_factor = mean_solution_arc_cost * intensification_lb;
    auto shaking_ub_factor = mean_solution_arc_cost * intensification_ub;

    const auto omega_base = std::max(1, static_cast<int>(std::ceil(std::log(instance.get_vertices_num()))));
    auto omega = std::vector<int>(instance.get_vertices_num(), omega_base);
    auto random_choice = std::uniform_int_distribution(0, 1);

#ifdef GUI
    auto renderer = Renderer(instance, best_solution.get_cost(), omega);
#endif


    const auto sa_initial_temperature = static_cast<float>(mean_arc_cost) * param.get_sa_initial_factor();
    const auto sa_final_temperature = static_cast<float>(mean_arc_cost) * param.get_sa_final_factor();

    auto sa = cobra::SimulatedAnnealing(sa_initial_temperature, sa_final_temperature, rand_engine, coreopt_iterations);

#ifdef VERBOSE
    std::cout << "Simulated annealing temperature goes from " << sa_initial_temperature << " to " << sa_final_temperature << ".\n\n";
#endif

    solution.commit();

    auto lk = LinKernighan(instance, tolerance);

    // auto pool = cobra::RouteSet(instance);

    for (auto iter = 0; iter < coreopt_iterations; iter++) {

        auto neighbor = solution;
        assert(neighbor.is_feasible());

#ifdef VERBOSE
        if (std::chrono::duration_cast<std::chrono::minutes>(std::chrono::high_resolution_clock::now() - global_time_begin).count() >= elapsed_minutes + 5) {
            printer.notify(
                "Optimizing for " +
                std::to_string(std::chrono::duration_cast<std::chrono::minutes>(std::chrono::high_resolution_clock::now() - global_time_begin).count()) +
                " minutes.");
            elapsed_minutes += 5;
        }
#endif


        auto removed = std::vector<int>();
        const auto walk_seed = rr.apply(neighbor, omega, removed);

        assert(neighbor.is_feasible());

#ifdef GUI
        const auto shaken_solution_cost = neighbor.get_cost();
#endif

        ruined_customers.clear();
        for (auto i = neighbor.get_cache_begin(); i != neighbor.get_cache_end(); i = neighbor.get_cache_next(i)) {
            ruined_customers.emplace_back(i);
        }

        cobra::randomly_flip_routes(instance, neighbor, rand_engine, 0.05f);

#ifdef VERBOSE
        welford_rac_after_shaking.update(static_cast<float>(neighbor.get_unstaged_changes().size()));
        welford_shaken_solutions.update(neighbor.get_cost());
#endif

        local_search.sequential_apply(neighbor);

        assert(neighbor.is_feasible());

#ifdef GUI
        const auto local_optimum_cost = neighbor.get_cost();
#endif

        average_number_of_vertices_accessed.update(static_cast<float>(neighbor.get_unstaged_changes().size()));
        const auto max_non_improving_iterations = static_cast<int>(
            std::ceil(delta * static_cast<float>(coreopt_iterations) * static_cast<float>(average_number_of_vertices_accessed.get_mean()) /
                      static_cast<float>(instance.get_vertices_num())));

#ifdef GUI
        if (iter % 100 == 0) {
            renderer.draw(best_solution, neighbor.get_unstaged_changes(), move_generators);
        }
#endif

#ifdef VERBOSE
        welford_rac_before_shaking.update(static_cast<float>(neighbor.get_unstaged_changes().size()));
        welford_local_optima.update(neighbor.get_cost());
#endif


        if (neighbor.get_cost() < best_solution.get_cost()) {


#ifdef VERBOSE
            printer.set_style(cobra::PrettyPrinter::BACKGROUND_CYAN);
            print_progress(printer, instance, iter, coreopt_iterations, main_opt_loop_begin_time, gamma, omega, sa, neighbor, param.get_bks_value()
    #ifdef LO_STORAGE
                                                                                                                                  ,
                           storage
    #endif
            );
            printer.unset_style();

            const auto old_cost = neighbor.get_cost();
#endif

            lk.apply(neighbor);
            assert(neighbor.is_feasible());

            // pool.add_routes(neighbor);

#ifdef VERBOSE
            const auto new_cost = neighbor.get_cost();

            if (new_cost + tolerance < old_cost) {

                printer.set_style(cobra::PrettyPrinter::BACKGROUND_BLUE);
                print_progress(printer, instance, iter, coreopt_iterations, main_opt_loop_begin_time, gamma, omega, sa, neighbor, param.get_bks_value()
    #ifdef LO_STORAGE
                                                                                                                                      ,
                               storage
    #endif
                );
                printer.unset_style();
            }
#endif

            best_solution = neighbor;

            gamma_vertices.clear();
            for (auto i = neighbor.get_cache_begin(); i != neighbor.get_cache_end(); i = neighbor.get_cache_next(i)) {
                gamma[i] = gamma_base;
                gamma_counter[i] = 0;
                gamma_vertices.emplace_back(i);
            }
            move_generators.set_active_percentage(gamma, gamma_vertices);

#ifdef VERBOSE
            welford_local_optima.reset();
            welford_local_optima.update(neighbor.get_cost());
            welford_shaken_solutions.reset();
            welford_shaken_solutions.update(neighbor.get_cost());
#endif

        } else {

            for (auto i = neighbor.get_cache_begin(); i != neighbor.get_cache_end(); i = neighbor.get_cache_next(i)) {
                gamma_counter[i]++;
                if (gamma_counter[i] >= max_non_improving_iterations) {
                    gamma[i] = std::min(gamma[i] * 2.0f, 1.0f);
                    gamma_counter[i] = 0;
                    gamma_vertices.clear();
                    gamma_vertices.emplace_back(i);
                    move_generators.set_active_percentage(gamma, gamma_vertices);
                }
            }
        }


        const auto seed_shake_value = omega[walk_seed];

        if (neighbor.get_cost() > shaking_ub_factor + solution.get_cost()) {
            for (auto i : ruined_customers) {
                if (omega[i] > seed_shake_value - 1) {
                    omega[i]--;
                }
            }

        } else if (neighbor.get_cost() >= solution.get_cost() && neighbor.get_cost() < solution.get_cost() + shaking_lb_factor) {
            for (auto i : ruined_customers) {
                if (omega[i] < seed_shake_value + 1) {
                    omega[i]++;
                }
            }

        } else {
            for (auto i : ruined_customers) {
                if (random_choice(rand_engine)) {
                    if (omega[i] > seed_shake_value - 1) {
                        omega[i]--;
                    }
                } else {
                    if (omega[i] < seed_shake_value + 1) {
                        omega[i]++;
                    }
                }
            }
        }

        assert(neighbor.is_feasible());
        if (sa.accept(solution, neighbor)) {

            solution = neighbor;

            assert(neighbor.is_feasible());
            assert(solution.is_feasible());

            if (!round_costs) {
                solution.recompute_costs();
            }  // avoid too many rounding errors get summed during LS

            solution.commit();
            const auto updated_mean_solution_arc_cost = solution.get_cost() / (static_cast<float>(instance.get_customers_num()) +
                                                                               2.0f * static_cast<float>(solution.get_routes_num()));
            shaking_lb_factor = updated_mean_solution_arc_cost * intensification_lb;
            shaking_ub_factor = updated_mean_solution_arc_cost * intensification_ub;
        }

        sa.decrease_temperature();

#ifdef GUI
        renderer.add_trajectory_point(shaken_solution_cost, local_optimum_cost, solution.get_cost(), solution.is_load_feasible(), best_solution.get_cost());
#endif

#ifdef VERBOSE
        partial_time_end = std::chrono::high_resolution_clock::now();
        const auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(partial_time_end - partial_time_begin).count();
        if (elapsed_time > 1) {
            print_progress(printer, instance, iter, coreopt_iterations, main_opt_loop_begin_time, gamma, omega, sa, best_solution, param.get_bks_value()
    #ifdef LO_STORAGE
                                                                                                                                       ,
                           storage
    #endif
            );
            partial_time_begin = std::chrono::high_resolution_clock::now();
        }
#endif
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