#ifndef _F4D_SPLITCHAIN_HPP_
#define _F4D_SPLITCHAIN_HPP_

#include "../SmallFlatMap.hpp"
#include "AbstractOperator.hpp"

namespace cobra {

    template <bool log_statistics>
    struct SplitChainStatistics { };
    template <>
    struct SplitChainStatistics<true> {

        Welford num_tree_nodes;
        Welford max_tree_depth;
    };

    template <bool handle_partial_solution = false, bool log_statistics = false, int max_split_nodes = 100>
    class SplitChain : SplitChainStatistics<log_statistics>, public AbstractOperator {

        static const int heap_unheaped = -1;
        static constexpr auto max_chain_length = 5;

        struct SplitNode {
            short heap_index = heap_unheaped;
            short predecessor = 0;
            short level = 0;
            bool direction = false;
            bool reversed = false;
            int infeas_demand = 0;
            int other_side_load = 0;
            float delta_sum = 0.0f;
            const MoveGenerator *move = nullptr;
            VerySmallFlatSet<int, 0, max_chain_length> visited_routes;
            // variables for resuming a paused search
            int i;
            int iNext;
            int next_movegen_index;
        };

        std::vector<SplitNode> split_nodes;
        std::vector<int> feasible_chains;

        BinaryHeapPtr<SplitNode, &SplitNode::heap_index, &SplitNode::infeas_demand> split_heap;
        inline int index_of(std::vector<SplitNode> &nodes, SplitNode &node) {
            return &node - nodes.data();
        }


    public:
        SplitChain(const Instance &instance_, MoveGenerators &moves_, float tolerance_) : AbstractOperator(instance_, moves_, tolerance_) {

            split_nodes.resize(max_split_nodes + 10);
            feasible_chains.reserve(max_split_nodes);
            // heap_array.resize(max_split_nodes + 10);
        }

        static constexpr bool is_symmetric = true;

    protected:
        inline void pre_processing(__attribute__((unused)) Solution &solution) override {

            for (int route = solution.get_first_route(); route != Solution::dummy_route; route = solution.get_next_route(route)) {
                solution.update_cumulative_route_loads(route);
            }
        }

        inline float compute_cost(const Solution &solution, const MoveGenerator &move) override {
            const auto i = move.get_first_vertex();
            const auto j = move.get_second_vertex();

            const auto iRoute = solution.get_route_index(i, j);
            const auto jRoute = solution.get_route_index(j, i);

            const auto iNext = solution.get_next_vertex(iRoute, i);
            const auto jNext = solution.get_next_vertex(jRoute, j);

            return -solution.get_next_cost(i) + this->instance.get_cost(i, j) - solution.get_next_cost(j) +
                   this->instance.get_cost(jNext, iNext);
        }

        bool is_feasible(const Solution &solution, const MoveGenerator &generating_move) override {

            int rni = 0;  // no. of generated nodes

            feasible_chains.clear();

            // First check whether the current `generating_move` is itself a feasible split move.
            // If this is the case, apply it without further searching for feasible ejection chains.

            auto capacity = this->instance.get_vehicle_capacity();
            auto min_capacity = this->instance.get_demand_sum() - solution.get_routes_num() * capacity;

            auto i = generating_move.get_first_vertex();
            auto j = generating_move.get_second_vertex();

            auto iRoute = solution.get_route_index(i, j);
            auto jRoute = solution.get_route_index(j, i);

            if (iRoute == jRoute) {
                return false;
            }

            assert(i != this->instance.get_depot());
            assert(j != this->instance.get_depot());

            // ensure iRoute feasibility
            const auto iRoute_load = solution.get_route_load_before_included(i) + solution.get_route_load_before_included(j);
            if (iRoute_load < min_capacity || iRoute_load > capacity) {
                return false;
            }

            // check jRoute feasibility
            const auto jRoute_load = solution.get_route_load_after_included(i) - this->instance.get_demand(i) + solution.get_route_load_after_included(j) -
                                     this->instance.get_demand(j);
            if (jRoute_load <= capacity) {

                // we found a feasible node! Setup some stuff used in the execute and return

                feasible_chains.push_back(0);
                split_nodes[0].move = &generating_move;
                split_nodes[0].predecessor = -1;
                split_nodes[rni].delta_sum = generating_move.get_delta();
                split_nodes[0].infeas_demand = 0;

                return true;
            }


            // do not start a chain if it is unlikely to close it
            const auto max_demand = 1.5 * capacity;
            if (jRoute_load >= max_demand) {
                return false;
            }

            // ... otherwise, start it!

            split_heap.reset();  // heap_reset();

            // we are now working on jRoute:
            // . iNext --- depot (originally in iRoute), we are going forward from iNext
            // . depot --- jNext (originally in jRoute), we are going backward from jNext and making sure to reverse get_next_vertex meaning because that route
            // part is now reversed
            //   Also in this node we may consider the the arc (jNext, iNext) which starts from jRoute and reaches iRoute

            const auto iNext = solution.get_next_vertex(i);

            if (iNext != this->instance.get_depot()) {  // is there a path to explore?
                split_nodes[rni].move = &generating_move;
                split_nodes[rni].delta_sum = generating_move.get_delta();
                split_nodes[rni].infeas_demand = jRoute_load - capacity;
                split_nodes[rni].level = 0;
                split_nodes[rni].visited_routes.clear();
                split_nodes[rni].visited_routes.insert(iRoute);
                split_nodes[rni].visited_routes.insert(jRoute);
                split_nodes[rni].predecessor = -1;
                split_nodes[rni].i = iNext;
                split_nodes[rni].iNext = solution.get_next_vertex(iNext);
                split_nodes[rni].direction = true;  // going forward
                split_nodes[rni].reversed = false;
                split_nodes[rni].other_side_load = solution.get_route_load_after_included(j) - this->instance.get_demand(j);
                split_nodes[rni].next_movegen_index = 0;
                split_heap.insert(&split_nodes[rni]);  // heap_insert(rni);
                rni++;
            }

            const auto jNext = solution.get_next_vertex(j);

            if (jNext != this->instance.get_depot()) {  // is there a path to explore?
                split_nodes[rni].move = &generating_move;
                split_nodes[rni].delta_sum = generating_move.get_delta();
                split_nodes[rni].infeas_demand = jRoute_load - capacity;
                split_nodes[rni].level = 0;
                split_nodes[rni].visited_routes.clear();
                split_nodes[rni].visited_routes.insert(iRoute);
                split_nodes[rni].visited_routes.insert(jRoute);
                split_nodes[rni].predecessor = -1;
                split_nodes[rni].i = jNext;
                split_nodes[rni].iNext = iNext;     // initially move to iRoute to check (jNext, iNext)
                split_nodes[rni].direction = true;  // going backward (which is forward in the original route)
                split_nodes[rni].reversed = true;
                split_nodes[rni].other_side_load = solution.get_route_load_after_included(i) - this->instance.get_demand(i);
                split_nodes[rni].next_movegen_index = 0;
                split_heap.insert(&split_nodes[rni]);  // heap_insert(rni);

                rni++;
            }

            while (!split_heap.empty()) {
                auto &curr = *split_heap.spy(0);
                const auto curr_index = index_of(split_nodes, curr);  // do not remove the vertex

                const auto done = curr.direction ? (curr.reversed ? search_split_exchange_from<true, true>(solution, curr_index, rni, max_demand)
                                                                  : search_split_exchange_from<true, false>(solution, curr_index, rni, max_demand))
                                                 : (curr.reversed ? search_split_exchange_from<false, true>(solution, curr_index, rni, max_demand)
                                                                  : search_split_exchange_from<false, false>(solution, curr_index, rni, max_demand));

                if (done) {
                    break;
                }
            }

        //end:

            if constexpr (log_statistics) {
                auto max_depth = 0;
                for (auto idx = 0; idx < rni; idx++) {
                    if (split_nodes[idx].level > max_depth) {
                        max_depth = split_nodes[idx].level;
                    }
                }
                this->max_tree_depth.update(max_depth);

                this->num_tree_nodes.update(rni);
            }

            return !feasible_chains.empty();
        }

        /**
         * The names used in this function are as if we are always going forward, just to make it simple to reason about.
         * The lambdas "transparently" manage the opposite direction.
         */
        template <bool direction, bool reversed>
        inline bool search_split_exchange_from(const Solution &solution, int curr_index, int &rni, int max_demand) {

            auto &curr = split_nodes[curr_index];

            auto capacity = this->instance.get_vehicle_capacity();
            auto min_capacity = this->instance.get_demand_sum() - solution.get_routes_num() * capacity;

            // solution.print(solution.get_route_index(curr.move->get_first_vertex()));
            // solution.print(solution.get_route_index(curr.move->get_second_vertex()));

            // Explore the current path
            // Note that we consider iRoute the route that we are exploring and jRoute the other one
            for (; curr.i != this->instance.get_depot();) {

                const auto i = curr.i;
                assert(i != this->instance.get_depot());

                // Keep the "visited" load updated
                const auto iDemand = this->instance.get_demand(i);
                curr.other_side_load += iDemand;

                const auto iNext = curr.iNext;  // we update this at the end

                const auto iDelta = this->instance.get_cost(i, iNext);

                while (curr.next_movegen_index < static_cast<int>(this->moves.get_move_generator_indices_involving_1st(i).size())) {

                    // make sure if we resume, we don't re-consider the same move gen
                    const auto move_index = this->moves.get_move_generator_indices_involving_1st(i)[curr.next_movegen_index++];

                    const auto &move = this->moves.get(move_index);
                    const auto j = move.get_second_vertex();

                    if constexpr (handle_partial_solution) {
                        if (!solution.is_vertex_in_solution(j)) {
                            continue;
                        }
                    }

                    if (j == this->instance.get_depot()) {
                        continue;
                    }

                    const auto jRoute = solution.get_route_index(j);

                    // skip forbidden routes
                    if (curr.visited_routes.count(jRoute) > 0) {
                        continue;
                    }

                    //const auto jDemand = this->instance.get_demand(j);

                    // Feasibility check for both routes

                    int iRoute_load;
                    int jRoute_load;

                    if constexpr (direction && reversed) {

                        // we are using solution.get_route_load_after_included(i) because the route is reversed
                        iRoute_load = solution.get_route_load_after_included(i) + solution.get_route_load_before_included(j);
                        jRoute_load = curr.other_side_load - this->instance.get_demand(i) + solution.get_route_load_after_included(j) -
                                      this->instance.get_demand(j);

                    } else if constexpr (direction && !reversed) {

                        iRoute_load = curr.other_side_load + solution.get_route_load_before_included(j);
                        const auto iNext_load = curr.iNext == this->instance.get_depot() ? 0 : solution.get_route_load_after_included(curr.iNext);
                        jRoute_load = iNext_load + solution.get_route_load_after_included(j) - this->instance.get_demand(j);

                    } else if constexpr (!direction && reversed) {

                        iRoute_load = curr.other_side_load + solution.get_route_load_before_included(j);

                        // iNext may be in another route!
                        const auto iNext_load = curr.iNext == this->instance.get_depot() ? 0 : solution.get_route_load_before_included(curr.iNext);
                        jRoute_load = iNext_load + solution.get_route_load_after_included(j) - this->instance.get_demand(j);

                    } else {  // !direction && !reversed

                        iRoute_load = solution.get_route_load_before_included(i) + solution.get_route_load_before_included(j);
                        jRoute_load = curr.other_side_load - this->instance.get_demand(i) + solution.get_route_load_after_included(j) -
                                      this->instance.get_demand(j);
                    }

                    const bool i_route_feasible = iRoute_load <= capacity;
                    const bool j_route_feasible = jRoute_load <= capacity;

                    const auto current_max_demand = std::max(iRoute_load, jRoute_load);
                    if (current_max_demand >= max_demand || std::min(iRoute_load, jRoute_load) < min_capacity) {
                        //((1.0 + 0.1 * max_chain_length) - 0.1 * std::min(max_chain_length, curr.level + 1)) * capacity) {
                        continue;
                    }
                    max_demand = current_max_demand;

                    // case 0: both infeasible
                    if (!(i_route_feasible || j_route_feasible)) {
                        continue;
                    }

                    // From now on we are sure to insert at least one SplitNode

                    const auto jNext = solution.get_next_vertex(j);  // jRoute is always traversed in the correct way
                    const auto delta = -iDelta + this->instance.get_cost(j, i) - solution.get_next_cost(j) + this->instance.get_cost(jNext, iNext);

                    // keep only improving chains
                    if (curr.delta_sum + delta > -this->tolerance) {
                        continue;
                    }


                    if (i_route_feasible && j_route_feasible) {
                        // boom! a feasible chain!

                        split_nodes[rni].move = &move;
                        split_nodes[rni].delta_sum = curr.delta_sum + delta;
                        split_nodes[rni].predecessor = curr_index;
                        split_nodes[rni].visited_routes = curr.visited_routes;
                        split_nodes[rni].visited_routes.insert(jRoute);
                        feasible_chains.push_back(rni);
                        return true;
                    }

                    if (curr.level + 1 >= max_chain_length - 1) {
                        continue;
                    }

                    if (rni >= max_split_nodes - 1) {
                        return true;
                    }

                    split_nodes[rni].next_movegen_index = 0;

                    if (j_route_feasible) {  // working on iRoute


                        if constexpr (direction && reversed) {

                            // i -> depot in iRoute with iNext = j in jRoute
                            split_nodes[rni].move = &move;
                            split_nodes[rni].delta_sum = curr.delta_sum + delta;
                            split_nodes[rni].infeas_demand = iRoute_load - capacity;
                            split_nodes[rni].level = curr.level + 1;
                            split_nodes[rni].predecessor = curr_index;
                            split_nodes[rni].visited_routes = curr.visited_routes;
                            split_nodes[rni].visited_routes.insert(jRoute);
                            split_nodes[rni].i = i;
                            split_nodes[rni].iNext = j;
                            split_nodes[rni].direction = true;
                            split_nodes[rni].reversed = true;
                            split_nodes[rni].other_side_load = solution.get_route_load_before_included(j);
                            split_nodes[rni].next_movegen_index = 0;
                            split_heap.insert(&split_nodes[rni]);  // heap_insert(rni);
                            rni++;
                            if (rni == max_split_nodes) {
                                return true;
                            }

                            // j -> depot in jRoute
                            split_nodes[rni].move = &move;
                            split_nodes[rni].delta_sum = curr.delta_sum + delta;
                            split_nodes[rni].infeas_demand = iRoute_load - capacity;
                            split_nodes[rni].level = curr.level + 1;
                            split_nodes[rni].predecessor = curr_index;
                            split_nodes[rni].visited_routes = curr.visited_routes;
                            split_nodes[rni].visited_routes.insert(jRoute);
                            split_nodes[rni].i = j;
                            split_nodes[rni].iNext = solution.get_prev_vertex(j);
                            split_nodes[rni].direction = false;
                            split_nodes[rni].reversed = true;
                            split_nodes[rni].other_side_load = solution.get_route_load_after_included(i);
                            split_nodes[rni].next_movegen_index = 0;
                            split_heap.insert(&split_nodes[rni]);  // heap_insert(rni);
                            rni++;
                            if (rni == max_split_nodes) {
                                return true;
                            }

                        } else if constexpr (direction && !reversed) {

                            // i -> j -> depot
                            split_nodes[rni].move = &move;
                            split_nodes[rni].delta_sum = curr.delta_sum + delta;
                            split_nodes[rni].infeas_demand = iRoute_load - capacity;
                            split_nodes[rni].level = curr.level + 1;
                            split_nodes[rni].predecessor = curr_index;
                            split_nodes[rni].visited_routes = curr.visited_routes;
                            split_nodes[rni].visited_routes.insert(jRoute);
                            split_nodes[rni].i = i;
                            split_nodes[rni].iNext = j;
                            split_nodes[rni].direction = false;
                            split_nodes[rni].reversed = true;
                            split_nodes[rni].other_side_load = curr.other_side_load - this->instance.get_demand(i);
                            split_nodes[rni].next_movegen_index = 0;
                            split_heap.insert(&split_nodes[rni]);  // heap_insert(rni);
                            rni++;
                            if (rni == max_split_nodes) {
                                return true;
                            }

                        } else if constexpr (!direction && reversed) {

                            // j -> depot
                            split_nodes[rni].move = &move;
                            split_nodes[rni].delta_sum = curr.delta_sum + delta;
                            split_nodes[rni].infeas_demand = iRoute_load - capacity;
                            split_nodes[rni].level = curr.level + 1;
                            split_nodes[rni].predecessor = curr_index;
                            split_nodes[rni].visited_routes = curr.visited_routes;
                            split_nodes[rni].visited_routes.insert(jRoute);
                            split_nodes[rni].i = i;
                            split_nodes[rni].iNext = j;
                            split_nodes[rni].direction = false;
                            split_nodes[rni].reversed = true;
                            split_nodes[rni].other_side_load = curr.other_side_load - this->instance.get_demand(i);
                            split_nodes[rni].next_movegen_index = 0;
                            split_heap.insert(&split_nodes[rni]);  // heap_insert(rni);
                            rni++;
                            if (rni == max_split_nodes) {
                                return true;
                            }

                        } else {  // !direction && !reversed

                            // i -> depot
                            split_nodes[rni].move = &move;
                            split_nodes[rni].delta_sum = curr.delta_sum + delta;
                            split_nodes[rni].infeas_demand = iRoute_load - capacity;
                            split_nodes[rni].level = curr.level + 1;
                            split_nodes[rni].predecessor = curr_index;
                            split_nodes[rni].visited_routes = curr.visited_routes;
                            split_nodes[rni].visited_routes.insert(jRoute);
                            split_nodes[rni].i = i;
                            split_nodes[rni].iNext = j;
                            split_nodes[rni].direction = false;
                            split_nodes[rni].reversed = false;
                            split_nodes[rni].other_side_load = solution.get_route_load_before_included(j);
                            split_nodes[rni].next_movegen_index = 0;
                            split_heap.insert(&split_nodes[rni]);  // heap_insert(rni);
                            rni++;
                            if (rni == max_split_nodes) {
                                return true;
                            }

                            // j -> depot
                            split_nodes[rni].move = &move;
                            split_nodes[rni].delta_sum = curr.delta_sum + delta;
                            split_nodes[rni].infeas_demand = iRoute_load - capacity;
                            split_nodes[rni].level = curr.level + 1;
                            split_nodes[rni].predecessor = curr_index;
                            split_nodes[rni].visited_routes = curr.visited_routes;
                            split_nodes[rni].visited_routes.insert(jRoute);
                            split_nodes[rni].i = j;
                            split_nodes[rni].iNext = solution.get_prev_vertex(j);
                            split_nodes[rni].direction = false;
                            split_nodes[rni].reversed = true;
                            split_nodes[rni].other_side_load = solution.get_route_load_before_included(i);
                            split_nodes[rni].next_movegen_index = 0;
                            split_heap.insert(&split_nodes[rni]);  // heap_insert(rni);
                            rni++;
                            if (rni == max_split_nodes) {
                                return true;
                            }
                        }

                    } else if (i_route_feasible) {  // working on jRoute

                        if constexpr (direction && reversed) {

                            if (jNext != this->instance.get_depot()) {

                                split_nodes[rni].move = &move;
                                split_nodes[rni].delta_sum = curr.delta_sum + delta;
                                split_nodes[rni].infeas_demand = jRoute_load - capacity;
                                split_nodes[rni].level = curr.level + 1;
                                split_nodes[rni].predecessor = curr_index;
                                split_nodes[rni].visited_routes = curr.visited_routes;
                                split_nodes[rni].visited_routes.insert(jRoute);
                                split_nodes[rni].i = jNext;
                                split_nodes[rni].iNext = iNext;
                                split_nodes[rni].direction = true;
                                split_nodes[rni].reversed = true;
                                split_nodes[rni].other_side_load = curr.other_side_load - this->instance.get_demand(i);
                                split_nodes[rni].next_movegen_index = 0;
                                split_heap.insert(&split_nodes[rni]);  // heap_insert(rni);
                                rni++;
                                if (rni == max_split_nodes) {
                                    return true;
                                }
                            }


                        } else if constexpr (direction && !reversed) {

                            if (iNext != this->instance.get_depot()) {

                                // iNext -> depot
                                split_nodes[rni].move = &move;
                                split_nodes[rni].delta_sum = curr.delta_sum + delta;
                                split_nodes[rni].infeas_demand = jRoute_load - capacity;
                                split_nodes[rni].level = curr.level + 1;
                                split_nodes[rni].predecessor = curr_index;
                                split_nodes[rni].visited_routes = curr.visited_routes;
                                split_nodes[rni].visited_routes.insert(jRoute);
                                split_nodes[rni].i = iNext;
                                split_nodes[rni].iNext = solution.get_next_vertex(iNext);
                                split_nodes[rni].direction = true;
                                split_nodes[rni].reversed = false;
                                split_nodes[rni].other_side_load = solution.get_route_load_after_included(j) - this->instance.get_demand(j);
                                split_nodes[rni].next_movegen_index = 0;
                                split_heap.insert(&split_nodes[rni]);  // heap_insert(rni);
                                rni++;
                                if (rni == max_split_nodes) {
                                    return true;
                                }
                            }

                            if (jNext != this->instance.get_depot()) {

                                // jNext -> depot
                                split_nodes[rni].move = &move;
                                split_nodes[rni].delta_sum = curr.delta_sum + delta;
                                split_nodes[rni].infeas_demand = jRoute_load - capacity;
                                split_nodes[rni].level = curr.level + 1;
                                split_nodes[rni].predecessor = curr_index;
                                split_nodes[rni].visited_routes = curr.visited_routes;
                                split_nodes[rni].visited_routes.insert(jRoute);
                                split_nodes[rni].i = jNext;
                                split_nodes[rni].iNext = iNext;
                                split_nodes[rni].direction = true;
                                split_nodes[rni].reversed = true;

                                // iNext may be in another route!
                                const auto iNext_load = curr.iNext == this->instance.get_depot() ? 0 : solution.get_route_load_after_included(curr.iNext);

                                split_nodes[rni].other_side_load = iNext_load;
                                split_nodes[rni].next_movegen_index = 0;
                                split_heap.insert(&split_nodes[rni]);  // heap_insert(rni);
                                rni++;
                                if (rni == max_split_nodes) {
                                    return true;
                                }
                            }


                        } else if constexpr (!direction && reversed) {


                            if (iNext != this->instance.get_depot()) {
                                split_nodes[rni].move = &move;
                                split_nodes[rni].delta_sum = curr.delta_sum + delta;
                                split_nodes[rni].infeas_demand = jRoute_load - capacity;
                                split_nodes[rni].level = curr.level + 1;
                                split_nodes[rni].predecessor = curr_index;
                                split_nodes[rni].visited_routes = curr.visited_routes;
                                split_nodes[rni].visited_routes.insert(jRoute);
                                split_nodes[rni].i = iNext;
                                split_nodes[rni].iNext = solution.get_prev_vertex(iNext);
                                split_nodes[rni].direction = false;
                                split_nodes[rni].reversed = true;
                                split_nodes[rni].other_side_load = solution.get_route_load_after_included(j) - this->instance.get_demand(j);
                                split_nodes[rni].next_movegen_index = 0;
                                split_heap.insert(&split_nodes[rni]);  // heap_insert(rni);

                                rni++;
                                if (rni == max_split_nodes) {
                                    return true;
                                }
                            }

                            if (jNext != this->instance.get_depot()) {
                                split_nodes[rni].move = &move;
                                split_nodes[rni].delta_sum = curr.delta_sum + delta;
                                split_nodes[rni].infeas_demand = jRoute_load - capacity;
                                split_nodes[rni].level = curr.level + 1;
                                split_nodes[rni].predecessor = curr_index;
                                split_nodes[rni].visited_routes = curr.visited_routes;
                                split_nodes[rni].visited_routes.insert(jRoute);
                                split_nodes[rni].i = jNext;
                                split_nodes[rni].iNext = iNext;
                                split_nodes[rni].direction = true;
                                split_nodes[rni].reversed = true;

                                // iNext may be in another route!
                                const auto iNext_load = curr.iNext == this->instance.get_depot() ? 0 : solution.get_route_load_before_included(curr.iNext);

                                split_nodes[rni].other_side_load = iNext_load;
                                split_nodes[rni].next_movegen_index = 0;
                                split_heap.insert(&split_nodes[rni]);  // heap_insert(rni);


                                rni++;
                                if (rni == max_split_nodes) {
                                    return true;
                                }
                            }

                        } else {  // !direction && !reversed

                            if (jNext != this->instance.get_depot()) {
                                split_nodes[rni].move = &move;
                                split_nodes[rni].delta_sum = curr.delta_sum + delta;
                                split_nodes[rni].infeas_demand = jRoute_load - capacity;
                                split_nodes[rni].level = curr.level + 1;
                                split_nodes[rni].predecessor = curr_index;
                                split_nodes[rni].visited_routes = curr.visited_routes;
                                split_nodes[rni].visited_routes.insert(jRoute);
                                split_nodes[rni].i = jNext;
                                split_nodes[rni].iNext = iNext;
                                split_nodes[rni].direction = true;
                                split_nodes[rni].reversed = true;
                                split_nodes[rni].other_side_load = curr.other_side_load - this->instance.get_demand(i);
                                split_nodes[rni].next_movegen_index = 0;
                                split_heap.insert(&split_nodes[rni]);  // heap_insert(rni);

                                rni++;
                                if (rni == max_split_nodes) {
                                    return true;
                                }
                            }
                        }
                    }

                    if (delta < -this->tolerance) {       // new nodes are more promising, get back to the heap and select them
                        curr.other_side_load -= iDemand;  // remove iDemand that will be re-added the node is re-retrieved
                        return false;
                    }
                }

                assert(curr.i != this->instance.get_depot());

                curr.next_movegen_index = 0;  // now we are changing vertex => consider all its move gen

                if constexpr (direction && reversed) {
                    curr.iNext = curr.i;
                    curr.i = solution.get_next_vertex(curr.i);
                    if (curr.i == this->instance.get_depot()) {
                        break;
                    }
                } else if constexpr (direction && !reversed) {
                    assert(curr.iNext == solution.get_next_vertex(curr.i));
                    curr.i = curr.iNext;
                    if (curr.i == this->instance.get_depot()) {
                        break;
                    }
                    curr.iNext = solution.get_next_vertex(curr.i);
                } else if constexpr (!direction && reversed) {
                    curr.i = curr.iNext;
                    if (curr.i == this->instance.get_depot()) {
                        break;
                    }
                    curr.iNext = solution.get_prev_vertex(curr.i);
                } else {  // !direction && !reversed
                    curr.iNext = curr.i;
                    curr.i = solution.get_prev_vertex(curr.i);
                    if (curr.i == this->instance.get_depot()) {
                        break;
                    }
                }
            }

            // remove curr from heap
            split_heap.get();  // heap_remove_top();

            return false;
        }

        inline void execute(Solution &solution, __attribute__((unused)) const MoveGenerator &p_move,
                            VertexSet &storage) override {

            assert(!feasible_chains.empty());

            // Apply the best chain after storing its vertices
            const auto best_chain_index = feasible_chains[0];

#ifndef NDEBUG
            auto &best_chain = split_nodes[best_chain_index];
            auto expected_delta = best_chain.delta_sum;
            auto old_cost = solution.get_cost();
#endif
            // TODO: manca la parte dove si riempie lo storage"

            auto split_moves = std::vector<int>();
            for (auto ptr = best_chain_index; ptr != -1; ptr = split_nodes[ptr].predecessor) {
                split_moves.emplace_back(ptr);
            }

            storage.insert(this->instance.get_depot());  // insert the depot once

            for (auto n = static_cast<int>(split_moves.size()) - 1; n >= 0; --n) {

                const auto ptr = split_moves[n];
                const auto move = split_nodes[ptr].move;

                const auto i = move->get_first_vertex();
                const auto j = move->get_second_vertex();

                // std::cout<< "Executing split from MG (" << i << ',' << j << ") Dir= " << split_nodes[ptr].direction << " Rev=" << split_nodes[ptr].reversed
                // << "\n";

                assert(i != this->instance.get_depot());
                assert(j != this->instance.get_depot());

                const auto iRoute = solution.get_route_index(i);
                const auto jRoute = solution.get_route_index(j);


                for (auto a = solution.get_first_customer(jRoute); a != this->instance.get_depot(); a = solution.get_next_vertex(a)) {
                    storage.insert(a);
                }
                storage.insert(i);

                // CUSTOM SPLIT (it may be optimized)
                const auto iNext = solution.get_next_vertex(i);
                //const auto jNext = solution.get_next_vertex(j);

                auto curr = j;
                while (curr != this->instance.get_depot()) {
                    const auto prev = solution.get_prev_vertex(curr);
                    solution.remove_vertex(jRoute, curr);
                    solution.insert_vertex_before(iRoute, iNext, curr);
                    curr = prev;
                }

                if (!solution.is_route_empty(jRoute)) {
                    solution.flip_route(jRoute);
                }

                curr = iNext;
                while (curr != this->instance.get_depot()) {
                    const auto next = solution.get_next_vertex(curr);
                    solution.remove_vertex(iRoute, curr);
                    solution.insert_vertex_before(jRoute, this->instance.get_depot(), curr);
                    curr = next;
                }


                if (solution.is_route_empty(iRoute)) {
                    solution.remove_route(iRoute);
                } else {
                    solution.update_cumulative_route_loads(iRoute);
                }

                if (solution.is_route_empty(jRoute)) {
                    solution.remove_route(jRoute);
                } else {
                    solution.update_cumulative_route_loads(jRoute);
                }
            }

            /*if (split_moves.size() > 1) {
                std::cout<< "Applied SplitChain of size " << moves.size() << "!!!!!\n";
                std::cout << "@@@ Remember we are checking solution feasibility after SPLIT @@@\n";
            }
            */

            if (!solution.is_feasible()) {
                abort();
            }

            assert(solution.is_feasible());

            assert(std::fabs(old_cost + expected_delta - solution.get_cost()) < 0.01f);
        }

        void post_processing(__attribute__((unused)) Solution &solution) override {
            // REMOVE
            // // std::cout<< " Called: " << called;
            // // std::cout<< ", Current FeasChainLenAVG: " << static_cast<double>(feasChainLenSum) / static_cast<double>(foudFeas);
            // // std::cout<< ", Current feasChainLenMax: " << feasChainLenMax;
            // // std::cout<< ", Current chainNumSumAVG: " << static_cast<double>(chainNumSum) / static_cast<double>(called);
            // // std::cout<< std::endl;
            // REMOVE
        }

        std::string get_additional_statistics() override {

            auto text = std::string();

            if constexpr (log_statistics) {
                text.append("MeanNumTreeNodes \t " + std::to_string(this->num_tree_nodes.get_mean()) + "\n");
                text.append("StdDevNumTreeNodes \t " + std::to_string(this->num_tree_nodes.get_standard_deviation()) + "\n");

                text.append("MeanNumTreeDepth \t " + std::to_string(this->max_tree_depth.get_mean()) + "\n");
                text.append("StdDevNumTreeDepth \t " + std::to_string(this->max_tree_depth.get_standard_deviation()) + "\n");
            }

            return text;
        }

        struct Cache12 {
            int v, next;
            float seqrem;
        };

        inline Cache12 prepare_cache12(const Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache12();
            c.v = vertex;
            c.next = solution.get_next_vertex(c.v);

            c.seqrem = -solution.get_next_cost(c.v);



            return c;
        }

        inline Cache12 prepare_cache12(const Solution &solution, int vertex, int backup) {

            auto c = Cache12();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v, backup);
            c.next = solution.get_next_vertex(route, c.v);

            c.seqrem = -solution.get_next_cost(c.v, c.next);



            return c;
        }


        inline float compute_cost(const struct Cache12 i, const struct Cache12 j) {

            const auto iSequenceAdd = this->instance.get_cost(i.v, j.v) + this->instance.get_cost(j.next, i.next);

            const auto delta = iSequenceAdd + i.seqrem + j.seqrem;

            return delta;
        }
    };

}  // namespace cobra

#endif