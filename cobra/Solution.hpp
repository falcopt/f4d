#ifndef _F4D_SOLUTION_HPP_
#define _F4D_SOLUTION_HPP_

#include <random>
#include <set>
#include <stack>
#include <unordered_set>
#include <vector>

#include "FixedSizeValueStack.hpp"
#include "Instance.hpp"
#include "LRUCache.hpp"
#include "macro.hpp"

namespace cobra {

    /**
     * Manage a CVRP solution.
     */
    class Solution {

        /**
         * Deep copy of a source solution.
         * @param source the source solution
         */
        void copy(const Solution &source) {

            routes_pool = source.routes_pool;

            depot_node = source.depot_node;
            for (int i = 0; i < instance.get_vertices_num(); i++) {
                customers_list[i] = source.customers_list[i];
            }
            for (int r = 0; r < instance.get_vertices_num(); r++) {
                routes_list[r] = source.routes_list[r];
            }
            solution_cost = source.solution_cost;

            cache = source.cache;
        }

    public:
        /**
         * Dummy value to identify whether a route contains real customers or it is empty.
         */
        static inline const int dummy_vertex = -1;

        /**
         * Dummy value to stop route's iteration. It is commonly used as
         * for(auto route = solution.get_first_route(); route != Solution::dummy_route; route = solution.get_next_route(route))
         */
        static inline const int dummy_route = 0;

        /**
         * Solution constructor. Note that the instance reference ownership is not transferred.
         * @param instance instance
         */
        explicit Solution(const Instance &instance) : Solution(instance, instance.get_vertices_num()) { }

        Solution(const Instance &instance, int history_len)
            : instance(instance),
              solution_cost(INFINITY),
              max_number_routes(instance.get_vertices_num() + 1),
              routes_pool(max_number_routes - 1,
                          [](int index) {
                              return index + 1;
                          }),
              depot_node({Solution::dummy_route, 0}),
              routes_list(new RouteNode[max_number_routes]),
              customers_list(new CustomerNode[instance.get_vertices_num()]),
              cache(history_len, instance.get_vertices_num()) { }

        /**
         * Generate a new solution from a the deep copy of a source solution.
         * @param source solution
         */
        Solution(const Solution &source)
            : instance(source.instance),
              solution_cost(INFINITY),
              max_number_routes(instance.get_vertices_num() + 1),
              routes_pool(max_number_routes - 1,
                          [](int index) {
                              return index + 1;
                          }),
              depot_node({Solution::dummy_route, 0}),
              routes_list(new RouteNode[max_number_routes]),
              customers_list(new CustomerNode[instance.get_vertices_num()]),
              cache(source.cache) {

            copy(source);
        }

        /**
         * Solution destructor.
         */
        virtual ~Solution() {
            delete[] customers_list;
            delete[] routes_list;
        }

        /**
         * Deep copy a source solution
         * @param source solution
         * @return solution
         */
        Solution &operator=(const Solution &source) {
            if (this == &source) {
                return *this;
            }
            copy(source);
            return *this;
        }

        /**
         * Reset a solution.
         */
        void reset() {

            solution_cost = 0.0;

            routes_pool.reset();

            depot_node.first_route = Solution::dummy_route;
            depot_node.num_routes = 0;

            for (int r = 0; r < max_number_routes; r++) {
                reset_route(r);
            }

            for (int i = 0; i < instance.get_vertices_num(); i++) {
                reset_vertex(i);
            }

            cache.clear();
        }

        /**
         * Returns the solution cost.
         * @return cost
         */
        float get_cost() const {
            return static_cast<float>(solution_cost);
        }

        /**
         * Returns the number of routes in the solution.
         * @return number of routes
         */
        inline auto get_routes_num() const -> int {
            return depot_node.num_routes;
        }

        /**
         * Build a one-customer route
         * @param customer
         * @return route index
         */
        int build_one_customer_route(const int customer) {

            assert(!is_customer_in_solution(customer));
            assert(customer != instance.get_depot());

            const auto route = request_route();

            customers_list[customer].prev = instance.get_depot();
            customers_list[customer].next = instance.get_depot();
            customers_list[customer].route_ptr = route;

            // head insert the route in the list
            const auto next_route = depot_node.first_route;
            routes_list[route].next = next_route;             // copy the previous head in the next
            depot_node.first_route = route;                   // and set the new head
            routes_list[route].prev = Solution::dummy_route;  // since it is an head insertion
            routes_list[next_route].prev = route;             // set the prev of the next route as this new route

            depot_node.num_routes++;

            routes_list[route].first_customer = customer;
            routes_list[route].last_customer = customer;
            routes_list[route].load = instance.get_demand(customer);
            routes_list[route].size = 1;

            solution_cost += 2.0f * instance.get_cost(instance.get_depot(), customer);

            cache.insert(customer);

            return route;
        }

        /**
         * Returns the index of the route serving a given customer.
         * @param customer
         * @return route index
         */
        inline int get_route_index(const int customer) const {
            assert(customer != instance.get_depot());
            return customers_list[customer].route_ptr;
        }

        /**
         * Returns the index if the route serving a given vertex. If vertex is the depot,
         * fallback is used to identify the route index.
         * @param vertex
         * @param fallback
         * @return route index
         */
        inline int get_route_index(const int vertex, const int fallback) const {
            if (unlikely(vertex == instance.get_depot())) {
                return customers_list[fallback].route_ptr;
            } else {
                return customers_list[vertex].route_ptr;
            }
        }

        /**
         * Returns the route load
         * @param route
         * @return load
         */
        inline int get_route_load(const int route) const {
            return routes_list[route].load;
        }

        /**
         * Index of the first route. Commonly used when iterating through all the routes
         * for(auto route = solution.get_first_route(); route != Solution::dummy_route; route = solution.get_next_route(route))
         * @return index of the first route
         */
        inline int get_first_route() const {
            return depot_node.first_route;
        }

        /**
         * Returns the next route. Commonly used when iterating through all the routes
         * for(auto route = solution.get_first_route(); route != Solution::dummy_route; route = solution.get_next_route(route))
         * @param route
         * @return index of the first route
         */
        inline int get_next_route(const int route) const {
            return routes_list[route].next;
        }

        inline int get_end_route() const {
            return dummy_route;
        }

        /**
         * Returns whether the route contains zero customers
         * @param route
         * @return boolean
         */
        inline bool is_route_empty(const int route) const {
            return routes_list[route].load == 0;
        }

        /**
         * Remove a vertex from a given route
         * @param route
         * @param vertex
         * @return the cost saving
         */
        float remove_vertex(const int route, const int vertex) {

            if (unlikely(vertex == instance.get_depot())) {

                const auto next = routes_list[route].first_customer;
                const auto prev = routes_list[route].last_customer;

                cache.insert(vertex);
                cache.insert(prev);
                cache.insert(next);

                set_prev_vertex_ptr(route, next, prev);
                set_next_vertex_ptr(route, prev, next);

                routes_list[route].first_customer = Solution::dummy_vertex;
                routes_list[route].last_customer = Solution::dummy_vertex;

                const auto delta = +instance.get_cost(prev, next) - instance.get_cost(prev, vertex) - instance.get_cost(vertex, next);

                solution_cost += delta;

                return delta;

            } else {

                assert(contains_vertex(route, vertex));
                assert(instance.get_depot() != vertex);

                const auto next = customers_list[vertex].next;
                const auto prev = customers_list[vertex].prev;

                cache.insert(vertex);
                cache.insert(prev);
                cache.insert(next);

                if (vertex == routes_list[route].first_customer) {
                    routes_list[route].first_customer = next;
                    set_prev_vertex_ptr(route, next, instance.get_depot());  // next might be the root of the route
                } else if (vertex == routes_list[route].last_customer) {
                    routes_list[route].last_customer = prev;
                    set_next_vertex_ptr(route, prev, instance.get_depot());  // prev might be the root of the route
                } else {
                    customers_list[prev].next = next;  // if vertex != route.first_customer => for sure prev is not the root
                    customers_list[next].prev = prev;  // if vertex != route.last_customer => for sure next is not the root
                }

                routes_list[route].load -= instance.get_demand(vertex);
                routes_list[route].size -= 1;

                const auto delta = +instance.get_cost(prev, next) - instance.get_cost(prev, vertex) - instance.get_cost(vertex, next);

                solution_cost += delta;

                // reset the removed vertex
                reset_vertex(vertex);

                return delta;
            }
        }

        /**
         * Removes an empty route. The route must be empty.
         * @param route empty route
         */
        void remove_route(const int route) {
            assert(is_route_empty(route));
            release_route(route);
        }

        /**
         * Returns the first customer of a route. It is commonly used as
         * for(auto i = solution.get_first_customer(route); i != instance.get_depot(); i = solution.get_next_vertex(i)){ ... }
         * @param route
         * @return first customer
         */
        inline auto get_first_customer(const int route) const -> int {
            return routes_list[route].first_customer;
        }

        /**
         * Returns the last customer of a route.
         * @param route
         * @return last customer
         */
        inline int get_last_customer(const int route) const {
            return routes_list[route].last_customer;
        }

        /**
         * Returns the vertex after a given one.
         * @param route
         * @param vertex
         * @return vertex after a given one
         */
        inline int get_next_vertex(const int route, const int vertex) const {

            assert(contains_vertex(route, vertex));

            if (unlikely(vertex == instance.get_depot())) {
                return routes_list[route].first_customer;
            } else {
                return customers_list[vertex].next;
            }
        }

        /**
         * Returns the vertex before a given one.
         * @param route
         * @param vertex
         * @return vertex before a given one
         */
        inline int get_prev_vertex(const int route, const int vertex) const {

            assert(contains_vertex(route, vertex));

            if (unlikely(vertex == instance.get_depot())) {
                return routes_list[route].last_customer;
            } else {
                return customers_list[vertex].prev;
            }
        }

        /**
         * Inserts a vertex before where in route.
         * @param route
         * @param where
         * @param vertex
         */
        void insert_vertex_before(const int route, const int where, const int vertex) {

            assert(where != vertex);

            if (unlikely(vertex == instance.get_depot())) {

                assert(routes_list[route].first_customer == Solution::dummy_vertex);
                assert(routes_list[route].last_customer == Solution::dummy_vertex);
                assert(where != instance.get_depot());

                const auto prev = customers_list[where].prev;

                cache.insert(prev);
                cache.insert(where);

                assert(prev != instance.get_depot());

                routes_list[route].first_customer = where;
                routes_list[route].last_customer = prev;

                customers_list[prev].next = instance.get_depot();
                customers_list[where].prev = instance.get_depot();

                const auto delta = +instance.get_cost(prev, instance.get_depot()) + instance.get_cost(instance.get_depot(), where) -
                                   instance.get_cost(prev, where);


                solution_cost += delta;

            } else {

                assert(!is_customer_in_solution(vertex));
                assert(vertex != instance.get_depot());

                const auto prev = get_prev_vertex(route, where);
                // insert vertex between prev and next

                cache.insert(prev);
                cache.insert(where);

                // vertex for sure is not the root of route
                customers_list[vertex].next = where;
                customers_list[vertex].prev = prev;
                customers_list[vertex].route_ptr = route;

                set_next_vertex_ptr(route, prev, vertex);
                set_prev_vertex_ptr(route, where, vertex);

                const auto delta = +instance.get_cost(prev, vertex) + instance.get_cost(vertex, where) - instance.get_cost(prev, where);

                solution_cost += delta;
                routes_list[route].load += instance.get_demand(vertex);
                routes_list[route].size += 1;
            }
        }

        /**
         * Reverses the sub-path identified by vertex_begin and vertex_end.
         * @param route
         * @param vertex_begin
         * @param vertex_end
         */
        void reverse_route_path(const int route, const int vertex_begin, const int vertex_end) {

            assert(vertex_begin != vertex_end);

            const auto pre = get_prev_vertex(route, vertex_begin);
            const auto stop = get_next_vertex(route, vertex_end);

            cache.insert(pre);
            cache.insert(stop);

            auto curr = vertex_begin;
            do {

                cache.insert(curr);

                const auto prev = get_prev_vertex(route, curr);
                const auto next = get_next_vertex(route, curr);

                set_prev_vertex_ptr(route, curr, next);
                set_next_vertex_ptr(route, curr, prev);

                curr = next;

            } while (curr != stop);

            if (vertex_end == pre && vertex_begin == stop) {
                // vertex_begin and vertex_end are contiguous
            } else {

                set_prev_vertex_ptr(route, vertex_end, pre);
                set_next_vertex_ptr(route, vertex_begin, stop);
                set_next_vertex_ptr(route, pre, vertex_end);
                set_prev_vertex_ptr(route, stop, vertex_begin);
            }

            const auto delta = -instance.get_cost(pre, vertex_begin) - instance.get_cost(vertex_end, stop) + instance.get_cost(pre, vertex_end) +
                               instance.get_cost(stop, vertex_begin);

            solution_cost += delta;
        }

        void flip_route(const int route) {

            const auto first = routes_list[route].first_customer;
            const auto last = routes_list[route].last_customer;

            routes_list[route].first_customer = last;
            routes_list[route].last_customer = first;

            auto curr = first;

            do {

                assert(curr != instance.get_depot());
                const auto prev = customers_list[curr].prev;
                const auto next = customers_list[curr].next;
                customers_list[curr].prev = next;
                customers_list[curr].next = prev;

                curr = next;

            } while (curr != instance.get_depot());

            assert(is_feasible(false));
        }

        /**
         * Merges two routes
         * @param route
         * @param route_to_append
         * @return index of the merge route
         */
        int append_route(const int route, const int route_to_append) {

            const auto route_end = routes_list[route].last_customer;
            const auto route_to_append_start = routes_list[route_to_append].first_customer;

            assert(route_end != instance.get_depot());
            assert(route_to_append_start != instance.get_depot());

            const auto delta = +instance.get_cost(route_end, route_to_append_start) - instance.get_cost(route_end, instance.get_depot()) -
                               instance.get_cost(instance.get_depot(), route_to_append_start);

            solution_cost += delta;

            customers_list[route_end].next = route_to_append_start;
            customers_list[route_to_append_start].prev = route_end;

            routes_list[route].last_customer = routes_list[route_to_append].last_customer;
            routes_list[route].load += routes_list[route_to_append].load;
            routes_list[route].size += routes_list[route_to_append].size;

            cache.insert(route_end);

            for (auto curr = route_to_append_start; curr != instance.get_depot(); curr = customers_list[curr].next) {
                customers_list[curr].route_ptr = route;

                cache.insert(curr);
            }

            release_route(route_to_append);

            return route;
        }

        /**
         * Generates a string representation of a given route.
         * @param route
         * @return string
         */
        std::string to_string(const int route) const {
            std::string str;
            str += "[" + std::to_string(route) + "] ";
            str += std::to_string(instance.get_depot()) + " ";
            for (int curr = routes_list[route].first_customer; curr != instance.get_depot(); curr = customers_list[curr].next) {
                str += std::to_string(curr) + " ";
            }
            str += std::to_string(instance.get_depot());
            return str;
        }

        /**
         * Prints a given route.
         * @param route
         */
        void print(const int route) const {
            if (is_missing_depot(route)) {
                std::cout << "Route " << route << " is in an INCONSISTENT state: missing the depot. It cannot be accessed without it.\n";
            } else {
                std::cout << to_string(route) << " (" << get_route_load(route) << ") " << get_route_cost(route) << "\n";
            }
        }

        /**
         * Prints the whole soluion.
         */
        void print() const {
            for (auto route = depot_node.first_route; route != Solution::dummy_route; route = routes_list[route].next) {
                print(route);
            }
            std::cout << "Solution cost = " << solution_cost << "\n";
        }

        void print_dimacs() const {

            for (auto route = get_first_route(), idx = 1; route != cobra::Solution::dummy_route; route = get_next_route(route), idx++) {
                std::cout << "Route #" << idx << ":";
                for (auto customer = get_first_customer(route); customer != instance.get_depot(); customer = get_next_vertex(customer)) {
                    std::cout << " " << customer;
                }
                std::cout << "\n";
            }

            std::cout << "Cost " << get_cost() << "\n";
            std::cout << std::flush;
        }

        /**
         * Updates the data structures containing the route's cumulative loads.
         * Commonly used in some local search operators (e.g. inter-route 2 opt)
         * @param route
         */
        void update_cumulative_route_loads(const int route) {

            assert(!is_route_empty(route));

            auto prev = routes_list[route].first_customer;

            customers_list[prev].load_before = instance.get_demand(prev);
            customers_list[prev].load_after = routes_list[route].load;

            auto curr = customers_list[prev].next;

            while (curr != instance.get_depot()) {

                customers_list[curr].load_before = customers_list[prev].load_before + instance.get_demand(curr);
                customers_list[curr].load_after = customers_list[prev].load_after - instance.get_demand(prev);

                prev = curr;
                curr = customers_list[curr].next;
            }
        }

        /**
         * Accesses the route's cumulative load before and including customer
         * @param customer
         * @return cumulative load
         */
        inline int get_route_load_before_included(const int customer) const {
            assert(customer != instance.get_depot());
            return customers_list[customer].load_before;
        }

        /**
         * Accesses the route's cumulative load after and including customer
         * @param customer
         * @return cumulative load
         */
        inline int get_route_load_after_included(const int customer) const {
            assert(customer != instance.get_depot());
            return customers_list[customer].load_after;
        }

        /**
         * Checks whether a route is currently used in the solution.
         * @param route
         * @return boolean flag
         */
        inline bool is_route_in_solution(const int route) const {
            return routes_list[route].first_customer != instance.get_depot() && routes_list[route].last_customer != instance.get_depot();
        }

        /**
         * Checks whether a customer is currently served in the solution.
         * @param customer
         * @return boolean flag
         */
        inline bool is_customer_in_solution(const int customer) const {
            assert(customer != instance.get_depot());
            return customers_list[customer].route_ptr != Solution::dummy_route;
        }

        inline bool is_vertex_in_solution(const int vertex) const {
            return vertex == instance.get_depot() || is_customer_in_solution(vertex);
        }

        /**
         * Check whether a vertex is served in a given route.
         * @param route
         * @param vertex
         * @return
         */
        inline bool contains_vertex(const int route, const int vertex) const {
            assert(vertex >= instance.get_vertices_begin() && vertex < instance.get_vertices_end() && route >= 0 && route < max_number_routes);
            return customers_list[vertex].route_ptr == route || vertex == instance.get_depot();
        }

        /**
         * Returns the number of customers served by the route
         * @param route
         * @return number of customers
         */
        inline int get_route_size(const int route) const {
            return routes_list[route].size;
        }

        /**
         * Returns the customer before a given one.
         * @param customer
         * @return customer before a given one
         */
        inline int get_prev_vertex(const int customer) const {
            assert(customer != instance.get_depot());
            return customers_list[customer].prev;
        }

        /**
         * Returns the customer after a given one. It is commonly used as
         * for(auto i = solution.get_first_customer(route); i != instance.get_depot(); i = solution.get_next_vertex(i)){ ... }
         * @param customer
         * @return customer after a given one
         */
        inline int get_next_vertex(const int customer) const {
            assert(customer != instance.get_depot());
            return customers_list[customer].next;
        }

        void swap_tails(const int i, const int iRoute, const int j, const int jRoute) {

            assert(i != instance.get_depot());
            assert(j != instance.get_depot());
            assert(iRoute != jRoute);
            assert(contains_vertex(iRoute, i));
            assert(contains_vertex(jRoute, j));

            /*const auto iNext = customers_list[i].next;

            auto curr = j;
            while (curr != instance.get_depot()) {
                const auto next = customers_list[curr].next;
                remove_vertex(jRoute, curr);
                insert_vertex_before(iRoute, iNext, curr);
                curr = next;
            }

            curr = iNext;
            while (curr != instance.get_depot()) {
                const auto next = customers_list[curr].next;
                remove_vertex(iRoute, curr);
                insert_vertex_before(jRoute, instance.get_depot(), curr);
                curr = next;
            }*/


            const auto jFirst = routes_list[jRoute].first_customer;
            const auto jLast = routes_list[jRoute].last_customer;
            const auto iLast = routes_list[iRoute].last_customer;

            const auto iNext = customers_list[i].next;
            const auto jPrev = customers_list[j].prev;


            // load update (we assume cumulative loads have been correctly updated)
            routes_list[iRoute].load = +customers_list[i].load_before + customers_list[j].load_after;
            routes_list[jRoute].load = customers_list[j].load_before - instance.get_demand(j) + customers_list[i].load_after - instance.get_demand(i);

            // cost update
            solution_cost += -instance.get_cost(i, iNext) - instance.get_cost(jPrev, j) + instance.get_cost(i, j) + instance.get_cost(jPrev, iNext);

            // update next, prev and route pointers
            customers_list[i].next = j;
            customers_list[j].prev = i;
            auto j_counter = 0;
            for (auto curr = j; curr != instance.get_depot(); curr = customers_list[curr].next) {
                customers_list[curr].route_ptr = iRoute;
                j_counter++;
            }

            set_next_vertex_ptr(jRoute, jPrev, iNext);
            set_prev_vertex_ptr(jRoute, iNext, jPrev);
            auto i_counter = 0;
            for (auto curr = iNext; curr != instance.get_depot(); curr = customers_list[curr].next) {
                customers_list[curr].route_ptr = jRoute;
                i_counter++;
            }

            // update routes size
            routes_list[iRoute].size += j_counter - i_counter;
            routes_list[jRoute].size += i_counter - j_counter;


            // update last and first customer if necessary

            // in the trivial case we just swap the tails
            routes_list[iRoute].last_customer = jLast;
            routes_list[jRoute].last_customer = iLast;

            if (j == jFirst) {
                routes_list[jRoute].first_customer = iNext;
                if (i == iLast) {
                    routes_list[jRoute].last_customer = instance.get_depot();
                };
            } else {
                if (i == iLast) {
                    routes_list[jRoute].last_customer = jPrev;
                }
            }

            cache.insert(i);
            cache.insert(j);
            cache.insert(iNext);
            cache.insert(jPrev);
        }

        void split(const int i, const int iRoute, const int j, const int jRoute) {

            assert(i != instance.get_depot());
            assert(j != instance.get_depot());

            const auto iNext = customers_list[i].next;
            const auto jNext = customers_list[j].next;

            auto curr = j;
            while (curr != instance.get_depot()) {
                const auto prev = customers_list[curr].prev;
                remove_vertex(jRoute, curr);
                insert_vertex_before(iRoute, iNext, curr);
                curr = prev;
            }

            auto before = jNext;
            curr = iNext;
            while (curr != instance.get_depot()) {
                const auto next = customers_list[curr].next;
                remove_vertex(iRoute, curr);
                insert_vertex_before(jRoute, before, curr);
                before = curr;
                curr = next;
            }
        }

        /**
         * Check whether the solution is CVRP feasible.
         * @return
         */
        bool is_feasible(const bool error_on_load_infeasible = true, const bool verbose = false) {

            std::vector<std::pair<std::string, int>> errors;
            std::vector<std::pair<std::string, int>> warnings;

            auto inconsistent_routes = std::set<int>();

            auto customers_already_visited_in_solution = std::set<int>();

            auto total_load = 0;
            auto total_cost = 0.0;

            // count how many times each vertex is a predecessor or successor of some other vertex
            auto predecessor_times = std::vector<int>(instance.get_vertices_num(), 0);
            auto successor_times = std::vector<int>(instance.get_vertices_num(), 0);

            for (auto route = get_first_route(); route != Solution::dummy_route; route = get_next_route(route)) {

                if (is_route_empty(route)) {
                    errors.emplace_back("Route " + std::to_string(route) + " is in solution but empty", __LINE__);
                }

                auto customers_already_visited_in_route = std::set<int>();

                auto initial_and_final_vertex = instance.get_depot();

                if (is_missing_depot(route)) {

                    warnings.emplace_back(
                        " Route " + std::to_string(route) +
                            " misses the depot. It is in an inconsistent state and there is no safe way to access it until depot is re-inserted",
                        __LINE__);
                    inconsistent_routes.insert(route);

                    // there is no simple way to access this route right now, the only way is to
                    // scan vertices to find a customer belonging to this route
                    for (auto c = instance.get_customers_begin(); c < instance.get_customers_end(); c++) {
                        if (customers_list[c].route_ptr == route) {
                            initial_and_final_vertex = c;
                            break;
                        }
                    }
                }

                auto route_load = 0;
                auto route_cost = 0.0;
                auto route_size = 0;

                auto curr = initial_and_final_vertex;
                auto next = Solution::dummy_vertex;
                do {

                    // check double visit
                    if (customers_already_visited_in_route.count(curr)) {
                        errors.emplace_back(
                            "Vertex " + std::to_string(curr) + " in route " + std::to_string(route) + " is visited more than once within this route", __LINE__);
                        break;
                    }

                    // check double visit
                    if (curr != instance.get_depot() && customers_already_visited_in_solution.count(curr)) {
                        errors.emplace_back(
                            "Vertex " + std::to_string(curr) + " in route " + std::to_string(route) + " is visited more than once in the solution", __LINE__);
                    }

                    next = get_next_vertex(route, curr);
                    const auto prev = get_prev_vertex(route, curr);

                    // check first customer
                    if (prev == instance.get_depot() && routes_list[route].first_customer != curr) {
                        errors.emplace_back("Vertex " + std::to_string(curr) + " in route " + std::to_string(route) +
                                                " has predecessor depot but it is not the first customer of the route which is instead vertex " +
                                                std::to_string(routes_list[route].first_customer),
                                            __LINE__);
                    }

                    // check last customer
                    if (next == instance.get_depot() && routes_list[route].last_customer != curr) {
                        errors.emplace_back("Vertex " + std::to_string(curr) + " in route " + std::to_string(route) +
                                                " has successor depot but it is not the last customer of the route which is instead vertex " +
                                                std::to_string(routes_list[route].last_customer),
                                            __LINE__);
                    }

                    // check linking pointers
                    if (curr != get_prev_vertex(route, next)) {
                        errors.emplace_back("Vertex " + std::to_string(curr) + " in route " + std::to_string(route) + " has successor " + std::to_string(next) +
                                                " but the predecessor of " + std::to_string(next) + " is instead vertex " +
                                                std::to_string(get_prev_vertex(route, next)),
                                            __LINE__);
                    }

                    if (curr != get_next_vertex(route, prev)) {
                        errors.emplace_back("Vertex " + std::to_string(curr) + " in route " + std::to_string(route) + " has predecessor " +
                                                std::to_string(prev) + " but the successor of " + std::to_string(prev) + " is instead vertex " +
                                                std::to_string(get_next_vertex(route, prev)),
                                            __LINE__);
                    }

                    if (curr != instance.get_depot() && customers_list[curr].route_ptr != route) {
                        errors.emplace_back("Vertex " + std::to_string(curr) + " in route " + std::to_string(route) + " has a route pointer " +
                                                std::to_string(customers_list[curr].route_ptr),
                                            __LINE__);
                    }

                    predecessor_times[prev]++;
                    successor_times[next]++;

                    customers_already_visited_in_route.insert(curr);
                    customers_already_visited_in_solution.insert(curr);

                    route_load += instance.get_demand(curr);
                    route_cost += instance.get_cost(curr, next);
                    if (curr != instance.get_depot()) {
                        route_size++;
                    }

                    curr = next;

                } while (next != initial_and_final_vertex);

                route_cost += instance.get_cost(curr, initial_and_final_vertex);  // return arc

                if (route_load != routes_list[route].load) {
                    errors.emplace_back("Route " + std::to_string(route) + " has a computed load of " + std::to_string(route_load) + " but the stored one is " +
                                            std::to_string(routes_list[route].load),
                                        __LINE__);
                }

                if (route_load > instance.get_vehicle_capacity()) {
                    const auto message = "Route " + std::to_string(route) + " has a load of " + std::to_string(route_load) + " but the vehicle capacity is " +
                                         std::to_string(instance.get_vehicle_capacity());
                    if (error_on_load_infeasible) {
                        errors.emplace_back(message, __LINE__);
                    } else {
                        warnings.emplace_back(message, __LINE__);
                    }
                }

                if (route_size != routes_list[route].size) {
                    errors.emplace_back("Route " + std::to_string(route) + " has a computed size of " + std::to_string(route_size) + " but the stored one is " +
                                            std::to_string(routes_list[route].size),
                                        __LINE__);
                }

                total_load += route_load;
                total_cost += route_cost;
            }

            if (predecessor_times[instance.get_depot()] != depot_node.num_routes) {
                errors.emplace_back("Depot is predecessor of " + std::to_string(predecessor_times[instance.get_depot()]) +
                                        " other vertices when it should be of exactly " + std::to_string(depot_node.num_routes),
                                    __LINE__);
            }

            if (successor_times[instance.get_depot()] != depot_node.num_routes) {
                errors.emplace_back("Depot is successor of " + std::to_string(successor_times[instance.get_depot()]) +
                                        " other vertices when it should be of exactly " + std::to_string(depot_node.num_routes),
                                    __LINE__);
            }

            for (auto i = instance.get_customers_begin(); i < instance.get_customers_end(); i++) {
                if (predecessor_times[i] > 1) {
                    errors.emplace_back("Vertex " + std::to_string(i) + " in route " + std::to_string(customers_list[i].route_ptr) + " is predecessor of " +
                                            std::to_string(predecessor_times[i]) + " other vertices when it should be of exactly 1",
                                        __LINE__);
                }
                if (successor_times[i] > 1) {
                    errors.emplace_back("Vertex " + std::to_string(i) + " in route " + std::to_string(customers_list[i].route_ptr) + " is successor of " +
                                            std::to_string(successor_times[i]) + " other vertices  when it should be of exactly 1",
                                        __LINE__);
                }
            }

            auto customers_not_served = std::vector<int>();
            auto customers_not_served_load = 0;
            for (auto i = instance.get_customers_begin(); i < instance.get_customers_end(); i++) {
                if (!is_customer_in_solution(i)) {
                    customers_not_served.push_back(i);
                    customers_not_served_load += instance.get_demand(i);
                } else {
                    if (routes_list[customers_list[i].route_ptr].first_customer <= instance.get_depot() ||
                        routes_list[customers_list[i].route_ptr].first_customer >= instance.get_vertices_end()) {
                        errors.emplace_back("Vertex " + std::to_string(i) + " in route " + std::to_string(customers_list[i].route_ptr) +
                                                " belongs to a route for which the first customer is " +
                                                std::to_string(routes_list[customers_list[i].route_ptr].first_customer),
                                            __LINE__);
                    }
                    if (routes_list[customers_list[i].route_ptr].last_customer <= instance.get_depot() ||
                        routes_list[customers_list[i].route_ptr].last_customer >= instance.get_vertices_end()) {
                        errors.emplace_back("Vertex " + std::to_string(i) + " in route " + std::to_string(customers_list[i].route_ptr) +
                                                " belongs to a route for which the last customer is " +
                                                std::to_string(routes_list[customers_list[i].route_ptr].last_customer),
                                            __LINE__);
                    }
                    const auto route = customers_list[i].route_ptr;
                    auto found = false;
                    for (auto curr = get_first_customer(route); curr != instance.get_depot(); curr = get_next_vertex(curr)) {
                        if (curr == i) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        errors.emplace_back("Vertex " + std::to_string(i) + " in route " + std::to_string(customers_list[i].route_ptr) +
                                                " cannot be found by scanning the route",
                                            __LINE__);
                    }
                }
            }

            if (!customers_not_served.empty()) {

                warnings.emplace_back("There are " + std::to_string(customers_not_served.size()) + " customers not served", __LINE__);
            }

            auto instance_demand_sum = instance.get_demand_sum();

            if (total_load + customers_not_served_load != instance_demand_sum) {
                errors.emplace_back("The load of served customers is " + std::to_string(total_load) + " but the sum of the load of all customers is " +
                                        std::to_string(instance_demand_sum),
                                    __LINE__);
            }

            if (std::fabs(get_cost() - total_cost) >= 0.5f) {
                errors.emplace_back(
                    "The solution has a computed cost of " + std::to_string(total_cost) + " but the stored one is " + std::to_string(get_cost()), __LINE__);
            }

            if (!errors.empty() || verbose) {
                std::cout << "== BEGIN OF SOLUTION FEASIBILITY CHECK REPORT ==\n";
            }

            if (!errors.empty()) {
                if (errors.size() == 1) {
                    std::cout << "There is 1 error\n";
                } else {
                    std::cout << "There are " << errors.size() << " errors\n";
                }
                for (const auto &entry : errors) {
                    const auto message = entry.first;
                    const auto line = entry.second;
                    std::cout << "+ LINE " << line << " + " << message << "\n";
                }
            }

            if (!errors.empty() || verbose) {
                if (warnings.size() == 1) {
                    std::cout << "There is 1 warning\n";
                } else {
                    std::cout << "There are " << warnings.size() << " warnings\n";
                }
                for (const auto &entry : warnings) {
                    const auto message = entry.first;
                    const auto line = entry.second;
                    std::cout << "+ LINE " << line << " + " << message << "\n";
                }
            }

            if (!errors.empty() || verbose) {
                std::cout << "== END OF SOLUTION FEASIBILITY CHECK REPORT ==\n";
            }

            return errors.empty();
        }

        /**
         * Returns the route cost.
         * @param route
         * @return route cost
         */
        float get_route_cost(const int route) const {
            auto curr = routes_list[route].first_customer;
            auto sum = instance.get_cost(instance.get_depot(), curr);
            while (curr != instance.get_depot()) {
                const auto next = customers_list[curr].next;
                sum += instance.get_cost(curr, next);
                curr = next;
            }

            return sum;
        }

        /**
         * Clear the set of recently modified vertices.
         */
        inline void commit() {
            cache.clear();
        }

        /**
         * Returns the recently modified vertices.
         * Note that it is not safe to perform any operations while iterating the returned set.
         * @return set
         */
        inline const LRUCache &get_unstaged_changes() const {
            return cache;
        }

        inline int get_cache_begin() const {
            return cache.begin();
        }

        inline int get_cache_next(const int i) const {
            return cache.get_next(i);
        }

        inline int get_cache_size() const {
            return cache.size();
        }

        inline int get_cache_end() const {
            return cache.end();
        }

        void set_unstaged_changes(LRUCache &cache) {
            cache = cache;
        }

        void recompute_costs() {

            solution_cost = 0.0;

            for (auto route = get_first_route(); route != dummy_route; route = get_next_route(route)) {

                solution_cost += instance.get_cost(instance.get_depot(), get_first_customer(route));
                auto curr = get_first_customer(route);
                do {
                    const auto next = get_next_vertex(curr);
                    solution_cost += instance.get_cost(curr, next);
                    curr = next;
                } while (curr != instance.get_depot());
            }
        }

        inline bool is_load_feasible(const int route) const {
            return routes_list[route].load <= instance.get_vehicle_capacity();
        }

        inline bool is_load_feasible() const {

            for (auto r = get_first_route(); r != dummy_route; r = get_next_route(r)) {
                if (!is_load_feasible(r)) {
                    return false;
                }
            }

            return true;
        }


    private:
        struct DepotNode {
            int first_route;
            int num_routes;
        };

        struct CustomerNode {
            int next;
            int prev;
            int route_ptr;
            int load_after;
            int load_before;
        };

        struct RouteNode {
            int first_customer;
            int last_customer;
            int load;
            int next;
            int prev;
            int size;
        };

        const cobra::Instance &instance;
        double solution_cost;
        const int max_number_routes;
        FixedSizeValueStack<int> routes_pool;
        struct DepotNode depot_node;
        struct RouteNode *routes_list;
        struct CustomerNode *customers_list;

        void reset_route(const int route) {
            routes_list[route].load = 0;
            routes_list[route].size = 0;
            routes_list[route].first_customer = instance.get_depot();
            routes_list[route].last_customer = instance.get_depot();
        }

        void reset_vertex(const int customer) {
            customers_list[customer].next = Solution::dummy_vertex;
            customers_list[customer].prev = Solution::dummy_vertex;
            customers_list[customer].route_ptr = Solution::dummy_route;
        }

        inline void set_next_vertex_ptr(const int route, const int vertex, const int next) {
            if (unlikely(vertex == instance.get_depot())) {
                routes_list[route].first_customer = next;
            } else {
                customers_list[vertex].next = next;
            }
        }

        inline void set_prev_vertex_ptr(const int route, const int vertex, const int prev) {
            if (unlikely(vertex == instance.get_depot())) {
                routes_list[route].last_customer = prev;
            } else {
                customers_list[vertex].prev = prev;
            }
        }

        inline int request_route() {

            assert(!routes_pool.is_empty());

            const auto route = routes_pool.get();

            return route;
        }

        inline void release_route(const int route) {

            const auto prevRoute = routes_list[route].prev;
            const auto nextRoute = routes_list[route].next;

            routes_list[prevRoute].next = nextRoute;
            routes_list[nextRoute].prev = prevRoute;
            depot_node.num_routes--;

            // head remove
            if (depot_node.first_route == route) {
                depot_node.first_route = nextRoute;
            }

            reset_route(route);

            routes_pool.push(route);
        }

        inline bool is_missing_depot(const int route) const {
            return get_first_customer(route) == Solution::dummy_vertex;
        }

        cobra::LRUCache cache;
    };


}  // namespace cobra


#endif