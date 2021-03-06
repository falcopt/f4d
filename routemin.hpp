#ifndef _F4D_ROUTEMIN_HPP_
#define _F4D_ROUTEMIN_HPP_

#include "cobra/LocalSearch.hpp"
#include "cobra/MoveGenerators.hpp"
#include "cobra/PrettyPrinter.hpp"
#include "cobra/Solution.hpp"


cobra::Solution routemin(const cobra::Instance &instance, const cobra::Solution &source, std::mt19937 &rand_engine, cobra::MoveGenerators &move_generators,
                         int kmin, int max_iter, float tolerance) {

    auto rvnd0 = cobra::RandomizedVariableNeighborhoodDescent<true>(
        instance, move_generators,
        {cobra::E11,   cobra::E10,  cobra::TAILS, cobra::SPLIT, cobra::RE22B, cobra::E22,  cobra::RE20,  cobra::RE21,  cobra::RE22S, cobra::E21, cobra::E20,
         cobra::TWOPT, cobra::RE30, cobra::E30,   cobra::RE33B, cobra::E33,   cobra::RE31, cobra::RE32B, cobra::RE33S, cobra::E31,   cobra::E32, cobra::RE32S},
        rand_engine, tolerance);

    auto local_search = cobra::VariableNeighborhoodDescentComposer(tolerance);
    local_search.append(&rvnd0);

    auto gamma_vertices = std::vector<int>();
    auto gamma = std::vector<float>(instance.get_vertices_num(), 1.0f);
    for (auto i = instance.get_vertices_begin(); i < instance.get_vertices_end(); i++) {
        gamma_vertices.emplace_back(i);
    }
    move_generators.set_active_percentage(gamma, gamma_vertices);

    auto best_solution = source;

    auto uniform_01_dist = std::uniform_real_distribution<float>(0.0f, 1.0f);
    auto customers_distribution = std::uniform_int_distribution(instance.get_customers_begin(), instance.get_customers_end() - 1);

    const auto t_base = 1.00f;
    const auto t_end = 0.01f;
    auto t = t_base;
    auto c = std::pow(t_end / t_base, 1.0 / max_iter);

    auto removed = std::vector<int>();
    removed.reserve(instance.get_customers_num());

    auto still_removed = std::vector<int>();
    still_removed.reserve(instance.get_customers_num());

    auto solution = best_solution;

    solution.commit();

    for (auto iter = 0; iter < max_iter; iter++) {

        // Remove all customers from the selected route and remove the route itself

        auto seed = cobra::Solution::dummy_vertex;
        do {
            seed = customers_distribution(rand_engine);
        } while (!solution.is_customer_in_solution(seed));
        auto selected_routes = std::vector<int>();
        selected_routes.push_back(solution.get_route_index(seed));
        const auto &neighbors = instance.get_neighbors_of(seed);


        for (auto n = 1u; n < neighbors.size(); n++) {
            const auto vertex = neighbors[n];
            if (vertex == instance.get_depot()) {
                continue;
            }
            if (!solution.is_customer_in_solution(vertex)) {
                continue;
            }
            const auto route = solution.get_route_index(vertex);
            if (route != selected_routes[0]) {
                selected_routes.push_back(route);
                break;
            }
        }

        removed.clear();
        removed.insert(removed.end(), still_removed.begin(), still_removed.end());
        still_removed.clear();

        for (auto selected_route : selected_routes) {

            auto curr = solution.get_first_customer(selected_route);
            do {
                const auto next = solution.get_next_vertex(curr);
                solution.remove_vertex(selected_route, curr);
                removed.emplace_back(curr);
                curr = next;
            } while (curr != instance.get_depot());

            solution.remove_route(selected_route);
        }

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

            for (auto route = solution.get_first_route(); route != cobra::Solution::dummy_route; route = solution.get_next_route(route)) {

                if (solution.get_route_load(route) + instance.get_demand(i) > instance.get_vehicle_capacity()) {
                    continue;
                }

                // before each customer
                for (auto j = solution.get_first_customer(route); j != instance.get_depot(); j = solution.get_next_vertex(j)) {
                    const auto prev = solution.get_prev_vertex(route, j);
                    const auto delta = -instance.get_cost(prev, j) + instance.get_cost(prev, i) + instance.get_cost(i, j);

                    if (delta < best_delta) {
                        best_route = route;
                        best_where = j;
                        best_delta = delta;
                    }
                }

                // before the depot
                const auto delta = -instance.get_cost(solution.get_last_customer(route), instance.get_depot()) +
                                   instance.get_cost(solution.get_last_customer(route), i) + instance.get_cost(i, instance.get_depot());

                if (delta < best_delta) {
                    best_route = route;
                    best_where = instance.get_depot();
                    best_delta = delta;
                }
            }

            if (best_route == -1) {

                const auto r = uniform_01_dist(rand_engine);

                if (r > t || solution.get_routes_num() < kmin) {
                    solution.build_one_customer_route(i);
                } else {
                    still_removed.push_back(i);
                }


            } else {
                solution.insert_vertex_before(best_route, best_where, i);
            }
        }

        local_search.sequential_apply(solution);

        solution.commit();

        if (still_removed.empty()) {

            if (solution.get_cost() < best_solution.get_cost() ||
                (solution.get_cost() == best_solution.get_cost() && solution.get_routes_num() < best_solution.get_routes_num())) {

                solution.print_dimacs();
                best_solution = solution;


                if (best_solution.get_routes_num() <= kmin) {
                    goto end;
                }
            }
        }

        const auto gap = (solution.get_cost() - best_solution.get_cost()) / best_solution.get_cost();

        if (gap > 0.0f) {

            solution = best_solution;
            still_removed.clear();
        }

        t *= c;

        assert(solution.is_feasible());
    }

end:

    assert(best_solution.is_feasible());

    return best_solution;
}

#endif