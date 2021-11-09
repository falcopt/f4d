#ifndef _F4D_COREOPTSOLVER_HPP_
#define _F4D_COREOPTSOLVER_HPP_

#include "LinKernighan.hpp"
#include "RuinAndRecreate.hpp"
#include "arg_parser.hpp"
#include "cobra/Instance.hpp"
#include "cobra/LocalSearch.hpp"
#include "cobra/MoveGenerators.hpp"
#include "cobra/NeighborAcceptance.hpp"
#include "cobra/PrettyPrinter.hpp"
#include "cobra/Solution.hpp"
#include "cobra/Welford.hpp"

#ifdef GUI
    #include "Renderer.hpp"
#endif

class CoreOptSolver {
public:
    CoreOptSolver(const cobra::Instance& instance_, const Parameters& params_, std::mt19937& rnd_, cobra::MoveGenerators& moves_,
                  cobra::VariableNeighborhoodDescentComposer& local_search_)
        : instance(instance_), param(params_), rand_engine(rnd_), rr(instance, rnd_), move_generators(moves_), local_search(local_search_) {

        mean_arc_cost = 0.0;
        for (auto i = instance.get_vertices_begin(); i < instance.get_vertices_end() - 1; i++) {
            for (auto j = i + 1; j < instance.get_vertices_end(); j++) {
                mean_arc_cost += instance.get_cost(i, j);
            }
        }
        mean_arc_cost /= (instance.get_vertices_num() * (instance.get_vertices_num() - 1) / 2.0);

        const auto omega_base = std::max(1, static_cast<int>(std::ceil(std::log(instance.get_vertices_num()))));
        omega = std::vector<int>(instance.get_vertices_num(), omega_base);
    }


    cobra::Solution fastopt(const cobra::Solution& warm_start, int coreopt_iterations) {


#ifdef VERBOSE
        auto partial_time_begin = std::chrono::high_resolution_clock::now();
        auto partial_time_end = std::chrono::high_resolution_clock::now();
#endif

        auto solution = warm_start;
        auto best_solution = solution;

        auto ruined_customers = std::vector<int>();

#ifdef VERBOSE
        std::cout << "Running FASTOPT for " << coreopt_iterations << " iterations.\n";

        auto printer = cobra::PrettyPrinter({
            {"%", cobra::PrettyPrinter::Field::Type::REAL, 5, " "},
            {"Iterations", cobra::PrettyPrinter::Field::Type::INTEGER, 10, " "},
            {param.get_bks_value().has_value() ? "Gap" : "Objective",
             param.needs_round_costs() && !param.get_bks_value().has_value() ? cobra::PrettyPrinter::Field::Type::INTEGER
                                                                             : cobra::PrettyPrinter::Field::Type::REAL,
             10, " ", 4},
            {"Routes", cobra::PrettyPrinter::Field::Type::INTEGER, 6, " "},
            {"Iter/s", cobra::PrettyPrinter::Field::Type::REAL, 10, " "},
            {"Eta (s)", cobra::PrettyPrinter::Field::Type::REAL, 10, " "},
            {"Omega", cobra::PrettyPrinter::Field::Type::REAL, 6, " "},
            {"Temp", cobra::PrettyPrinter::Field::Type::REAL, 6, " "},
        });

        auto main_opt_loop_begin_time = std::chrono::high_resolution_clock::now();

#endif

        const auto intensification_lb = param.get_shaking_lb_factor();
        const auto intensification_ub = param.get_shaking_ub_factor();

        const auto mean_solution_arc_cost = solution.get_cost() /
                                            (static_cast<float>(instance.get_customers_num()) + 2.0f * static_cast<float>(solution.get_routes_num()));

        auto shaking_lb_factor = mean_solution_arc_cost * intensification_lb;
        auto shaking_ub_factor = mean_solution_arc_cost * intensification_ub;

        auto random_choice = std::uniform_int_distribution(0, 1);

        const auto sa_initial_temperature = static_cast<float>(mean_arc_cost) * param.get_sa_initial_factor();
        const auto sa_final_temperature = static_cast<float>(mean_arc_cost) * param.get_sa_final_factor();

        auto sa = cobra::SimulatedAnnealing(sa_initial_temperature, sa_final_temperature, rand_engine, coreopt_iterations);

#ifdef GUI
        auto renderer = Renderer(instance, best_solution.get_cost(), omega);
#endif


        for (auto iter = 0; iter < coreopt_iterations; iter++) {

            auto neighbor = solution;
            assert(neighbor.is_feasible());

            auto removed = std::vector<int>();
            const auto walk_seed = rr.apply(neighbor, omega, removed);

            assert(neighbor.is_feasible());

#ifdef GUI
            const auto shaken_solution_cost = neighbor.get_cost();
            const auto local_optimum_cost = neighbor.get_cost();
#endif

            ruined_customers.clear();
            for (auto i = neighbor.get_cache_begin(); i != neighbor.get_cache_end(); i = neighbor.get_cache_next(i)) {
                ruined_customers.emplace_back(i);
            }

#ifdef GUI
            if (iter % 1000 == 0) {
                renderer.draw(best_solution, neighbor.get_unstaged_changes(), move_generators);
            }
#endif

            if (neighbor.get_cost() < best_solution.get_cost()) {
                best_solution = neighbor;
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

                if (!param.needs_round_costs()) {
                    solution.recompute_costs();
                }  // avoid too many rounding errors get summed up

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
                print_progress_fast(printer, iter, coreopt_iterations, main_opt_loop_begin_time, sa, best_solution, param.get_bks_value());
                partial_time_begin = std::chrono::high_resolution_clock::now();
            }
#endif
        }

        return best_solution;
    }

    cobra::Solution coreopt(const cobra::Solution& warm_start, int coreopt_iterations) {

#ifdef VERBOSE
        const auto global_time_begin = std::chrono::high_resolution_clock::now();
        auto partial_time_begin = std::chrono::high_resolution_clock::now();
        auto partial_time_end = std::chrono::high_resolution_clock::now();
#endif

        auto solution = warm_start;
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

        auto printer = cobra::PrettyPrinter({
            {"%", cobra::PrettyPrinter::Field::Type::REAL, 5, " "},
            {"Iterations", cobra::PrettyPrinter::Field::Type::INTEGER, 10, " "},
            {param.get_bks_value().has_value() ? "Gap" : "Objective",
             param.needs_round_costs() && !param.get_bks_value().has_value() ? cobra::PrettyPrinter::Field::Type::INTEGER
                                                                             : cobra::PrettyPrinter::Field::Type::REAL,
             10, " ", 4},
            {"Routes", cobra::PrettyPrinter::Field::Type::INTEGER, 6, " "},
            {"Iter/s", cobra::PrettyPrinter::Field::Type::REAL, 10, " "},
            {"Eta (s)", cobra::PrettyPrinter::Field::Type::REAL, 10, " "},
            {"Gamma", cobra::PrettyPrinter::Field::Type::REAL, 5, " "},
            {"Omega", cobra::PrettyPrinter::Field::Type::REAL, 6, " "},
            {"Temp", cobra::PrettyPrinter::Field::Type::REAL, 6, " "},
        });

        auto main_opt_loop_begin_time = std::chrono::high_resolution_clock::now();

        auto elapsed_minutes = 0;
#endif

        const auto intensification_lb = param.get_shaking_lb_factor();
        const auto intensification_ub = param.get_shaking_ub_factor();

        const auto mean_solution_arc_cost = solution.get_cost() /
                                            (static_cast<float>(instance.get_customers_num()) + 2.0f * static_cast<float>(solution.get_routes_num()));

        auto shaking_lb_factor = mean_solution_arc_cost * intensification_lb;
        auto shaking_ub_factor = mean_solution_arc_cost * intensification_ub;

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

        auto lk = LinKernighan(instance, param.get_tolerance());

        // auto pool = cobra::RouteSet(instance);

        for (auto iter = 0; iter < coreopt_iterations; iter++) {

            auto neighbor = solution;
            assert(neighbor.is_feasible());

#ifdef VERBOSE
            if (std::chrono::duration_cast<std::chrono::minutes>(std::chrono::high_resolution_clock::now() - global_time_begin).count() >=
                elapsed_minutes + 5) {
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


            if (neighbor.get_cost() < best_solution.get_cost()) {


#ifdef VERBOSE
                printer.set_style(cobra::PrettyPrinter::BACKGROUND_CYAN);
                print_progress(printer, iter, coreopt_iterations, main_opt_loop_begin_time, gamma, sa, neighbor, param.get_bks_value());
                printer.unset_style();
                const auto old_cost = neighbor.get_cost();
#endif

                lk.apply(neighbor);
                assert(neighbor.is_feasible());

                // pool.add_routes(neighbor);

#ifdef VERBOSE
                const auto new_cost = neighbor.get_cost();

                if (new_cost + param.get_tolerance() < old_cost) {

                    printer.set_style(cobra::PrettyPrinter::BACKGROUND_BLUE);
                    print_progress(printer, iter, coreopt_iterations, main_opt_loop_begin_time, gamma, sa, neighbor, param.get_bks_value());
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

                if (!param.needs_round_costs()) {
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
                print_progress(printer, iter, coreopt_iterations, main_opt_loop_begin_time, gamma, sa, best_solution, param.get_bks_value());
                partial_time_begin = std::chrono::high_resolution_clock::now();
            }
#endif
        }

        return best_solution;
    }

private:
    const cobra::Instance& instance;
    const Parameters& param;
    std::mt19937& rand_engine;
    RuinAndRecreate rr;
    cobra::MoveGenerators& move_generators;
    cobra::VariableNeighborhoodDescentComposer& local_search;
    double mean_arc_cost;
    std::vector<int> omega;

    void print_progress_fast(cobra::PrettyPrinter& printer, int iter, const int iterations,
                             std::chrono::time_point<std::chrono::high_resolution_clock>& main_opt_loop_begin_time, cobra::SimulatedAnnealing& sa,
                             cobra::Solution& best_solution, std::optional<float> bks) {


        const auto progress = 100.0f * static_cast<float>(iter + 1) / static_cast<float>(iterations);
        const auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - main_opt_loop_begin_time)
                                         .count();
        const auto iter_per_second = static_cast<float>(iter + 1) / (static_cast<float>(elapsed_seconds) + 0.01f);
        const auto remaining_iter = iterations - iter;
        const auto estimated_rem_time = static_cast<float>(remaining_iter) / iter_per_second;

        auto omega_mean = 0.0f;
        for (auto i = instance.get_customers_begin(); i < instance.get_customers_end(); i++) {
            omega_mean += omega[i];
        }
        omega_mean /= static_cast<float>(instance.get_customers_num());


        printer.print(progress, iter + 1, bks.has_value() ? 100.0f * (best_solution.get_cost() - *bks) / *bks : best_solution.get_cost(),
                      best_solution.get_routes_num(), iter_per_second, estimated_rem_time, omega_mean, sa.get_temperature());
    }

    void print_progress(cobra::PrettyPrinter& printer, int iter, const int iterations,
                        std::chrono::time_point<std::chrono::high_resolution_clock>& main_opt_loop_begin_time, std::vector<float>& gamma,
                        cobra::SimulatedAnnealing& sa, cobra::Solution& best_solution, std::optional<float> bks) {


        const auto progress = 100.0f * static_cast<float>(iter + 1) / static_cast<float>(iterations);
        const auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - main_opt_loop_begin_time)
                                         .count();
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
};

#endif