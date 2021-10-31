#ifndef _F4D_ROUTEMIN_HPP_
#define _F4D_ROUTEMIN_HPP_

#include "cobra/LocalSearch.hpp"
#include "cobra/MoveGenerators.hpp"
#include "cobra/PrettyPrinter.hpp"
#include "cobra/Solution.hpp"
#include "cobra/NeighborAcceptance.hpp"

#ifdef GUI
    #include "Renderer.hpp"
#endif

cobra::Solution routemin(const cobra::Instance &instance, const cobra::Solution &source, std::mt19937 &rand_engine, cobra::MoveGenerators &move_generators,
                         int kmin, int max_iter, float tolerance, bool round_costs) {

#ifdef VERBOSE
    auto partial_time_begin = std::chrono::high_resolution_clock::now();
    auto partial_time_end = std::chrono::high_resolution_clock::now();
#endif

    auto rvnd0 = cobra::RandomizedVariableNeighborhoodDescent<true, false>(
        instance, move_generators,
        {cobra::E11,   cobra::E10,  cobra::TAILS, cobra::SPLIT, cobra::RE22B, cobra::E22,  cobra::RE20,  cobra::RE21,  cobra::RE22S, cobra::E21, cobra::E20,
         cobra::TWOPT, cobra::RE30, cobra::E30,   cobra::RE33B, cobra::E33,   cobra::RE31, cobra::RE32B, cobra::RE33S, cobra::E31,   cobra::E32, cobra::RE32S},
        rand_engine, tolerance);

    auto local_search = cobra::VariableNeighborhoodDescentComposer(tolerance);
    local_search.append(&rvnd0);


    auto rvnd1 = cobra::RandomizedVariableNeighborhoodDescent(
        instance, move_generators,
        {cobra::E11,   cobra::E10,  cobra::TAILS, cobra::SPLIT, cobra::RE22B, cobra::E22,  cobra::RE20,  cobra::RE21,  cobra::RE22S, cobra::E21, cobra::E20,
         cobra::TWOPT, cobra::RE30, cobra::E30,   cobra::RE33B, cobra::E33,   cobra::RE31, cobra::RE32B, cobra::RE33S, cobra::E31,   cobra::E32, cobra::RE32S},
        rand_engine, tolerance);
    auto rvnd2 = cobra::RandomizedVariableNeighborhoodDescent(instance, move_generators, {cobra::EJCH, cobra::TLCH, cobra::STCH}, rand_engine, tolerance);

    auto local_search_feas = cobra::VariableNeighborhoodDescentComposer(tolerance);
    local_search_feas.append(&rvnd1);
    local_search_feas.append(&rvnd2);

    auto gamma_vertices = std::vector<int>();
    auto gamma = std::vector<float>(instance.get_vertices_num(), 0.25f);
    for (auto i = instance.get_vertices_begin(); i < instance.get_vertices_end(); i++) {
        gamma_vertices.emplace_back(i);
    }
    move_generators.set_active_percentage(gamma, gamma_vertices);

    auto best_solution = source;

    auto uniform_01_dist = std::uniform_real_distribution<float>(0.0f, 1.0f);
    auto customers_distribution = std::uniform_int_distribution(instance.get_customers_begin(), instance.get_customers_end() - 1);

    const auto t_base = 1.0f;
    const auto t_end = 0.01f;
    auto t = t_base;
    auto c = std::pow(t_end / t_base, 1.0 / max_iter);

    auto removed = std::vector<int>();
    removed.reserve(instance.get_customers_num());

    auto solution_still_removed = std::vector<int>();
    solution_still_removed.reserve(instance.get_customers_num());


    auto solution = best_solution;

    solution.commit();


#ifdef VERBOSE
    const auto main_opt_loop_begin_time = std::chrono::high_resolution_clock::now();

    auto printer = cobra::PrettyPrinter({{"%", cobra::PrettyPrinter::Field::Type::INTEGER, 3, " "},
                                         {"Objective", cobra::PrettyPrinter::Field::Type::INTEGER, 10, " "},
                                         {"Routes", cobra::PrettyPrinter::Field::Type::INTEGER, 6, " "},
                                         {"Iter/s", cobra::PrettyPrinter::Field::Type::REAL, 7, " "},
                                         {"Eta (s)", cobra::PrettyPrinter::Field::Type::REAL, 6, " "},
                                         {"% Inf", cobra::PrettyPrinter::Field::Type::REAL, 6, " "}});

    auto number_infeasible_solutions = 0;

#endif


#ifdef GUI

    auto served = std::vector<int>(instance.get_vertices_num());
    for (auto i = instance.get_vertices_begin(); i < instance.get_vertices_end(); i++) {
        served[i] = 1;
    }

    auto renderer = Renderer(instance, best_solution.get_cost(), served);
#endif

    auto mean_arc_cost = 0.0;
    for (auto i = instance.get_vertices_begin(); i < instance.get_vertices_end() - 1; i++) {
        for (auto j = i + 1; j < instance.get_vertices_end(); j++) {
            mean_arc_cost += instance.get_cost(i, j);
        }
    }
    mean_arc_cost /= (instance.get_vertices_num() * (instance.get_vertices_num() - 1) / 2.0);


    const auto sa_initial_temperature = mean_arc_cost * 100;  // 10.0f
    const auto sa_final_temperature = mean_arc_cost / 10.0f;  // 100.0f

    auto sa = cobra::SimulatedAnnealing(sa_initial_temperature, sa_final_temperature, rand_engine, max_iter);

    for (auto iter = 0; iter < max_iter; iter++) {

#ifdef VERBOSE
        partial_time_end = std::chrono::high_resolution_clock::now();
        const auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(partial_time_end - partial_time_begin).count();
        if (elapsed_time > 1) {

            const auto progress = 100.0f * (iter + 1.0f) / static_cast<float>(max_iter);
            const auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - main_opt_loop_begin_time)
                                             .count();
            const auto iter_per_second = static_cast<float>(iter + 1.0f) / (static_cast<float>(elapsed_seconds) + 0.01f);
            const auto remaining_iter = max_iter - iter;
            const auto estimated_rem_time = static_cast<float>(remaining_iter) / iter_per_second;
            const auto fraction_infeasible_solutions = static_cast<float>(number_infeasible_solutions) / (iter + 1.0f);

            printer.print(progress, best_solution.get_cost(), best_solution.get_routes_num(), iter_per_second, estimated_rem_time,
                          fraction_infeasible_solutions);


            partial_time_begin = std::chrono::high_resolution_clock::now();
        }
#endif

        auto neighbor = solution;
        auto neighbor_still_removed = solution_still_removed;

        // Remove all customers from the selected route and remove the route itself

        auto seed = cobra::Solution::dummy_vertex;
        do {
            seed = customers_distribution(rand_engine);
        } while (!neighbor.is_customer_in_solution(seed));
        auto selected_routes = std::vector<int>();
        selected_routes.push_back(neighbor.get_route_index(seed));
        const auto &neighbors = instance.get_neighbors_of(seed);


        for (auto n = 1u; n < neighbors.size(); n++) {
            const auto vertex = neighbors[n];
            if (vertex == instance.get_depot()) {
                continue;
            }
            if (!neighbor.is_customer_in_solution(vertex)) {
                continue;
            }
            const auto route = neighbor.get_route_index(vertex);
            if (route != selected_routes[0]) {
                selected_routes.push_back(route);
                break;
            }
        }


        removed.clear();
        removed.insert(removed.end(), neighbor_still_removed.begin(), neighbor_still_removed.end());
        neighbor_still_removed.clear();

        for (auto selected_route : selected_routes) {

            auto curr = neighbor.get_first_customer(selected_route);
            do {
                const auto next = neighbor.get_next_vertex(curr);
                neighbor.remove_vertex(selected_route, curr);
                removed.emplace_back(curr);
                curr = next;
            } while (curr != instance.get_depot());

            neighbor.remove_route(selected_route);
        }

#ifdef GUI
        const auto shaken_solution_cost = neighbor.get_cost();
#endif

        if (rand_engine() % 2 == 0) {
            std::sort(removed.begin(), removed.end(), [&instance](auto i, auto j) {
                return instance.get_demand(i) > instance.get_demand(j);
            });
        } else {
            std::shuffle(removed.begin(), removed.end(), rand_engine);
        }


        for (auto i : removed) {

            auto best_route = -1;
            auto best_where = -1;
            auto best_delta = std::numeric_limits<float>::max();

            for (auto route = neighbor.get_first_route(); route != neighbor.get_end_route(); route = neighbor.get_next_route(route)) {

                if (neighbor.get_route_load(route) + instance.get_demand(i) > instance.get_vehicle_capacity()) {
                    continue;
                }

                // before each customer
                for (auto j = neighbor.get_first_customer(route); j != instance.get_depot(); j = neighbor.get_next_vertex(j)) {
                    const auto prev = neighbor.get_prev_vertex(route, j);
                    const auto delta = -instance.get_cost(prev, j) + instance.get_cost(prev, i) + instance.get_cost(i, j);

                    if (delta < best_delta
#ifdef DISTANCE_CONSTRAINED
                        && neighbor.get_route_cost(route) + delta <= instance.get_max_route_cost()
#endif
                    ) {
                        best_route = route;
                        best_where = j;
                        best_delta = delta;
                    }
                }

                // before the depot
                const auto delta = -instance.get_cost(neighbor.get_last_customer(route), instance.get_depot()) +
                                   instance.get_cost(neighbor.get_last_customer(route), i) + instance.get_cost(i, instance.get_depot());

                if (delta < best_delta
#ifdef DISTANCE_CONSTRAINED
                    && neighbor.get_route_cost(route) + delta <= instance.get_max_route_cost()
#endif
                ) {
                    best_route = route;
                    best_where = instance.get_depot();
                    best_delta = delta;
                }
            }

            if (best_route == -1) {

                const auto r = uniform_01_dist(rand_engine);

                if (r > t || neighbor.get_routes_num() < kmin - 1) {
                    neighbor.build_one_customer_route(i);
                } else {
                    neighbor_still_removed.push_back(i);
                }


            } else {

                neighbor.insert_vertex_before(best_route, best_where, i);
            }
        }

        if (neighbor_still_removed.empty()) {
            local_search_feas.sequential_apply(neighbor);
        } else {
            local_search.sequential_apply(neighbor);
        }


#ifdef GUI
        const auto local_optimum_cost = neighbor.get_cost();
#endif


#ifdef GUI
        if (iter % 100 == 0) {
            for (auto i = instance.get_vertices_begin(); i < instance.get_vertices_end(); i++) {
                served[i] = static_cast<int>(!neighbor.is_vertex_in_solution(i));
            }
            renderer.draw(best_solution, neighbor.get_unstaged_changes(), move_generators);
        }
#endif

        solution.commit();

        if (neighbor_still_removed.empty()) {

            if (neighbor.get_cost() < best_solution.get_cost() ||
                (neighbor.get_cost() == best_solution.get_cost() && neighbor.get_routes_num() < best_solution.get_routes_num())) {
                // if(neighbor.get_routes_num() < best_solution.get_routes_num() || (neighbor.get_routes_num() == best_solution.get_routes_num() &&
                // neighbor.get_cost() < best_solution.get_cost())) {
                best_solution = neighbor;

#ifdef DIMACS
                best_solution.print_dimacs();
#endif

                if (best_solution.get_routes_num() <= kmin) {
                    goto end;
                }
            }

        } else {

#ifdef VERBOSE
            number_infeasible_solutions++;
#endif
        }

#ifdef GUI
        renderer.add_trajectory_point(shaken_solution_cost, local_optimum_cost, solution.get_cost(), neighbor_still_removed.empty(), best_solution.get_cost());
#endif

        if (sa.accept(solution, neighbor)) {

            solution = neighbor;
            solution_still_removed = neighbor_still_removed;

            solution.commit();

            if (!round_costs) {
                solution.recompute_costs();
            }  // avoid too many rounding errors get summed during LS

            assert(neighbor.is_feasible());
            assert(solution.is_feasible());
        }

        t *= c;
        sa.decrease_temperature();

        assert(solution.is_feasible());
    }

end:

    assert(best_solution.is_feasible());

    return best_solution;
}

#endif