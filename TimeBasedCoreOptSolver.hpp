#ifndef _F4D_TIMEBASEDCOREOPTSOLVER_HPP_
#define _F4D_TIMEBASEDCOREOPTSOLVER_HPP_

#include "LinKernighan.hpp"
#include "RuinAndRecreate.hpp"
#include "arg_parser.hpp"
#include "cobra/Instance.hpp"
#include "cobra/LocalSearch.hpp"
#include "cobra/MoveGenerators.hpp"
#include "cobra/NeighborAcceptance.hpp"
#include "cobra/PrettyPrinter.hpp"
#include "cobra/SPH.hpp"
#include "cobra/Solution.hpp"
#include "cobra/Welford.hpp"

class TimeBasedCoreOptSolver {
public:
    TimeBasedCoreOptSolver(const cobra::Instance& instance_, const Parameters& params_, std::mt19937& rnd_, cobra::MoveGenerators& moves_,
                           cobra::VariableNeighborhoodDescentComposer& local_search_, sph::SPHeuristic& sph_)
        : TimeBasedCoreOptSolver(instance_, params_, rnd_, moves_, local_search_) {
        sph = std::addressof(sph_);
    }

    TimeBasedCoreOptSolver(const cobra::Instance& instance_, const Parameters& params_, std::mt19937& rnd_, cobra::MoveGenerators& moves_,
                           cobra::VariableNeighborhoodDescentComposer& local_search_)
        : instance(instance_), param(params_), rand_engine(rnd_), rr(instance, rnd_), move_generators(moves_), local_search(local_search_), sph(nullptr) {

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

    cobra::Solution coreopt(const cobra::Solution& warm_start, int coreopt_iterations) {

        const auto global_time_begin = std::chrono::high_resolution_clock::now();

        auto solution = warm_start;
        auto best_solution = solution;

        std::vector<sph::real_t> routes_cost;

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

        const auto intensification_lb = param.get_shaking_lb_factor();
        const auto intensification_ub = param.get_shaking_ub_factor();

        const auto mean_solution_arc_cost = solution.get_cost() /
                                            (static_cast<float>(instance.get_customers_num()) + 2.0f * static_cast<float>(solution.get_routes_num()));

        auto shaking_lb_factor = mean_solution_arc_cost * intensification_lb;
        auto shaking_ub_factor = mean_solution_arc_cost * intensification_ub;

        auto random_choice = std::uniform_int_distribution(0, 1);

        const auto sa_initial_temperature = static_cast<float>(mean_arc_cost) * param.get_sa_initial_factor();
        const auto sa_final_temperature = static_cast<float>(mean_arc_cost) * param.get_sa_final_factor();

        auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - global_time_begin).count();
        auto sa = cobra::TimeBasedSimulatedAnnealing(sa_initial_temperature, sa_final_temperature, rand_engine, coreopt_iterations - elapsed_time);


        solution.commit();

        auto lk = LinKernighan(instance, param.get_tolerance());

        for (auto iter = 0; elapsed_time < coreopt_iterations; iter++) {

            auto neighbor = solution;
            assert(neighbor.is_feasible());

            auto removed = std::vector<int>();
            const auto walk_seed = rr.apply(neighbor, omega, removed);

            ruined_customers.clear();
            for (auto i = neighbor.get_cache_begin(); i != neighbor.get_cache_end(); i = neighbor.get_cache_next(i)) {
                ruined_customers.emplace_back(i);
            }

            cobra::randomly_flip_routes(instance, neighbor, rand_engine, 0.05f);

            local_search.sequential_apply(neighbor);

            average_number_of_vertices_accessed.update(static_cast<float>(neighbor.get_unstaged_changes().size()));

            const auto iter_per_second = static_cast<float>(iter + 1) / (static_cast<float>(elapsed_time) + 0.01f);
            const auto remaining_time = coreopt_iterations - elapsed_time;
            const auto estimated_remaining_iter = iter_per_second * remaining_time;
            const auto expected_total_iterations_num = iter + 1 + estimated_remaining_iter;

            const auto max_non_improving_iterations = static_cast<int>(
                std::ceil(delta * static_cast<float>(expected_total_iterations_num) * static_cast<float>(average_number_of_vertices_accessed.get_mean()) /
                          static_cast<float>(instance.get_vertices_num())));


            routes_cost.clear();
            for (auto route = neighbor.get_first_route(); route != cobra::Solution::dummy_route; route = neighbor.get_next_route(route))
                routes_cost.push_back(neighbor.get_route_cost(route));
            sph->add_solution(routes_cost, neighbor.get_routes());


            if (neighbor.get_cost() < best_solution.get_cost()) {

                lk.apply(neighbor);
                neighbor.print_dimacs();

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

            elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - global_time_begin).count();
            if (sa.accept(solution, neighbor, elapsed_time)) {

                solution = neighbor;

                if (!param.needs_round_costs()) {
                    solution.recompute_costs();
                }  // avoid too many rounding errors get summed during LS

                solution.commit();
                const auto updated_mean_solution_arc_cost = solution.get_cost() / (static_cast<float>(instance.get_customers_num()) +
                                                                                   2.0f * static_cast<float>(solution.get_routes_num()));
                shaking_lb_factor = updated_mean_solution_arc_cost * intensification_lb;
                shaking_ub_factor = updated_mean_solution_arc_cost * intensification_ub;
            }
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
    sph::SPHeuristic* sph;
};

#endif