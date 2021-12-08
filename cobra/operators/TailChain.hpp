#ifndef _F4D_TAILCHAIN_HPP_
#define _F4D_TAILCHAIN_HPP_

#include <unordered_set>

#include "../SmallFlatMap.hpp"
#include "AbstractOperator.hpp"


namespace cobra {

    template <bool log_statistics>
    struct TailChainStatistics { };
    template <>
    struct TailChainStatistics<true> {

        Welford num_tree_nodes;
        Welford max_tree_depth;
    };

    template <bool handle_partial_solution = false, bool log_statistics = false, int max_tail_nodes = 100>
    class TailChain : TailChainStatistics<log_statistics>, public AbstractOperator {

        static const int heap_unheaped = -1;
        static constexpr auto max_chain_length = 5;

        // REMOVE
        // unsigned called = 0;
        // unsigned foudFeas = 0;
        // long long unsigned feasChainLenSum = 0;
        // unsigned feasChainLenMax = 0;
        // double averageLengthSum = 0.0;
        // long long unsigned chainNumSum = 0;
        // REMOVE

        struct TailNode {
            short heap_index = heap_unheaped;
            short predecessor = 0;
            short level = 0;
            bool direction = false;
            int infeas_demand = 0;
            int other_side_load = 0;
            float delta_sum = 0.0f;
            const MoveGenerator *move = nullptr;
            VerySmallFlatSet<int, 0, max_chain_length> visited_routes;

            int i;
            int iNext;
            int next_movegen_index;
        };

        std::vector<TailNode> tail_nodes;
        int feasible_chain;

        /* HEAP SECTION */

        BinaryHeapPtr<TailNode, &TailNode::heap_index, &TailNode::infeas_demand> tail_heap;
        inline int index_of(std::vector<TailNode> &nodes, TailNode &node) {
            return &node - nodes.data();
        }

    public:
        TailChain(const Instance &instance_, MoveGenerators &moves_, float tolerance_) : AbstractOperator(instance_, moves_, tolerance_) {

            tail_nodes.resize(max_tail_nodes + 10);
            feasible_chain = -1;
            // heap_array.resize(max_tail_nodes + 10);
        }

        static constexpr bool is_symmetric = false;

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
            const auto jPrev = solution.get_prev_vertex(jRoute, j);

            return -solution.get_next_cost(i) + this->instance.get_cost(i, j) - solution.get_prev_cost(j) +
                   this->instance.get_cost(jPrev, iNext);
        }

        bool is_feasible(const Solution &solution, const MoveGenerator &generating_move) override {

            int rni = 0;  // no. of generated nodes

            feasible_chain = -1;

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
            int i_route_load = solution.get_route_load_before_included(i) + solution.get_route_load_after_included(j);
            if (i_route_load < min_capacity || i_route_load > capacity) {
                return false;
            }

            // check for jRoute feasibility
            int j_route_load = solution.get_route_load_before_included(j) - this->instance.get_demand(j) + solution.get_route_load_after_included(i) -
                               this->instance.get_demand(i);
            if (j_route_load <= capacity) {
                feasible_chain = 0;
                tail_nodes[0].move = &generating_move;
                tail_nodes[0].delta_sum = generating_move.get_delta();
                tail_nodes[0].infeas_demand = 0;
                tail_nodes[0].predecessor = -1;
                return true;
            }

            const auto max_demand = 1.5 * capacity;
            if (j_route_load >= max_demand) {
                return false;
            }
            // If the generating move is not feasible by itself, we start a chain

            tail_heap.reset();  // heap_reset();

            // set up state variables

            auto iNext = solution.get_next_vertex(iRoute, i);
            if (iNext != this->instance.get_depot()) {
                tail_nodes[rni].move = &generating_move;
                tail_nodes[rni].delta_sum = generating_move.get_delta();
                tail_nodes[rni].infeas_demand = j_route_load - capacity;
                tail_nodes[rni].level = 0;
                tail_nodes[rni].visited_routes.clear();
                tail_nodes[rni].visited_routes.insert(iRoute);
                tail_nodes[rni].visited_routes.insert(jRoute);
                tail_nodes[rni].predecessor = -1;
                tail_nodes[rni].i = iNext;
                tail_nodes[rni].iNext = solution.get_next_vertex(iNext);
                tail_nodes[rni].direction = true;
                tail_nodes[rni].other_side_load = solution.get_route_load_before_included(j) - this->instance.get_demand(j);
                tail_nodes[rni].next_movegen_index = 0;
                tail_heap.insert(&tail_nodes[rni]);  // heap_insert(rni);

                rni++;
            }

            auto jPrev = solution.get_prev_vertex(j);

            if (const auto jPrevPrev = solution.get_prev_vertex(jRoute, jPrev);
                jPrev != this->instance.get_depot() && jPrevPrev != this->instance.get_depot()) {
                tail_nodes[rni].move = &generating_move;
                tail_nodes[rni].delta_sum = generating_move.get_delta();
                tail_nodes[rni].infeas_demand = j_route_load - capacity;
                tail_nodes[rni].level = 0;
                tail_nodes[rni].visited_routes.clear();
                tail_nodes[rni].visited_routes.insert(iRoute);
                tail_nodes[rni].visited_routes.insert(jRoute);
                tail_nodes[rni].predecessor = -1;
                tail_nodes[rni].i = jPrev;
                tail_nodes[rni].iNext = iNext;
                tail_nodes[rni].direction = false;
                tail_nodes[rni].other_side_load = solution.get_route_load_after_included(i) - this->instance.get_demand(i);
                tail_nodes[rni].next_movegen_index = 0;
                tail_heap.insert(&tail_nodes[rni]);  // heap_insert(rni);

                rni++;
            }

            while (!tail_heap.empty()) {
                auto &curr = *tail_heap.spy(0);
                const auto curr_index = index_of(tail_nodes, curr);  // do not remove the vertex

                const auto done = curr.direction ? search_tail_exchange_from<true>(solution, curr_index, rni, max_demand)
                                                 : search_tail_exchange_from<false>(solution, curr_index, rni, max_demand);

                if (done) {
                    break;
                }
            }

            // end:
            if constexpr (log_statistics) {
                auto max_depth = 0;
                for (auto idx = 0; idx < rni; idx++) {
                    if (tail_nodes[idx].level > max_depth) {
                        max_depth = tail_nodes[idx].level;
                    }
                }
                this->max_tree_depth.update(max_depth);

                this->num_tree_nodes.update(rni);
            }

            // REMOVE
            //++called;
            // int count = 0;
            // if (feasible_chain >= 0) {
            //    ++foudFeas;
            //    for (auto ptr = feasible_chain; ptr != -1; ptr = relocation_nodes[ptr].predecessor)
            //        ++count;
            //}
            // feasChainLenSum += count;
            // if (feasChainLenMax < count)
            //    feasChainLenMax = count;
            // chainNumSum += rni;
            // REMOVE
            // std::cout << " return: end of is_feasible()\n";
            return feasible_chain >= 0;
        }

        /**
         * The names used in this function are as if we are always going forward, just to make it simple to reason about.
         * The lambdas "transparently" manage the opposite direction.
         */
        template <bool direction>
        inline bool search_tail_exchange_from(const Solution &solution, int curr_index, int &rni, int max_demand) {

            auto capacity = this->instance.get_vehicle_capacity();
            auto min_capacity = this->instance.get_demand_sum() - solution.get_routes_num() * capacity;

            // Convenience Lambda to generalize for both directions
            /*const auto next_vertex = [&](TailNode &curr) {
                if constexpr (direction) {
                    return curr.iNext;
                } else {
                    return solution.get_prev_vertex(curr.i);
                }
            };*/

            auto &curr = tail_nodes[curr_index];

            while (curr.i != this->instance.get_depot()) {

                const auto i = curr.i;
                assert(i != this->instance.get_depot());

                // Keep the "visited" load updated
                const auto iDemand = this->instance.get_demand(i);
                curr.other_side_load += iDemand;

                const auto iNext = curr.iNext;

                // So: we need the load after iNext because we cannot use i since it might be of the wrong route. However iNext can be the depot so, let's
                // hanlde the case here and stop thinking about it. The important thing is that we never compute the load after i directly.
                const auto load_after_iNext = iNext == this->instance.get_depot() ? 0 : solution.get_route_load_after_included(iNext);
                const auto iPrev = solution.get_prev_vertex(i);
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

                    assert(j != this->instance.get_depot());
                    const auto jRoute = solution.get_route_index(j);

                    // skip forbidden routes
                    if (curr.visited_routes.count(jRoute) > 0) {
                        continue;
                    }

                    const auto jDemand = this->instance.get_demand(j);

                    // Feasibility check for both routes

                    const auto [i_route_load, j_route_load] = [&] {
                        if constexpr (direction) {
                            return std::make_pair(curr.other_side_load + solution.get_route_load_after_included(j),
                                                  load_after_iNext + solution.get_route_load_before_included(j) - jDemand);
                        } else {
                            return std::make_pair(solution.get_route_load_before_included(i) + solution.get_route_load_after_included(j),
                                                  curr.other_side_load + solution.get_route_load_before_included(j) - iDemand - jDemand);
                        }
                    }();

                    const bool i_route_feasible = i_route_load <= capacity;
                    const bool j_route_feasible = j_route_load <= capacity;

                    const auto current_max_demand = std::max(i_route_load, j_route_load);
                    if (current_max_demand >= max_demand || std::min(i_route_load, j_route_load) < min_capacity) {
                        //((1.0 + 0.1 * max_chain_length) - 0.1 * std::min(max_chain_length, curr.level + 1)) * capacity) {
                        continue;
                    }
                    max_demand = current_max_demand;

                    // case 0: both infeasible
                    if (!(i_route_feasible || j_route_feasible)) {
                        continue;
                    }

                    // From now on we are sure to insert at least one TailNode
                    const auto jPrev = solution.get_prev_vertex(j);
                    const auto delta = -iDelta + this->instance.get_cost(i, j) - solution.get_prev_cost(j) + this->instance.get_cost(jPrev, iNext);

                    // keep only improving chains
                    if (curr.delta_sum + delta > -this->tolerance) {
                        continue;
                    }

                    // common fields
                    tail_nodes[rni].move = &move;
                    tail_nodes[rni].delta_sum = curr.delta_sum + delta;
                    tail_nodes[rni].level = curr.level + 1;
                    tail_nodes[rni].predecessor = curr_index;
                    tail_nodes[rni].visited_routes = curr.visited_routes;
                    tail_nodes[rni].visited_routes.insert(jRoute);
                    tail_nodes[rni].next_movegen_index = 0;


                    if (i_route_feasible && j_route_feasible) {  // boom! a feasible chain!
                        feasible_chain = rni;
                        return true;
                    }

                    if (curr.level + 1 >= max_chain_length - 1) {
                        continue;
                    }  // This chain cannot be completed, so skip it

                    if (rni >= max_tail_nodes - 1) {
                        return true;
                    }  // reached hard limit in the number of explored nodes

                    if (j_route_feasible) {
                        tail_nodes[rni].i = i;
                        tail_nodes[rni].iNext = j;
                        tail_nodes[rni].direction = true;
                        tail_nodes[rni].infeas_demand = i_route_load - capacity;


                        if constexpr (direction) tail_nodes[rni].other_side_load = curr.other_side_load - iDemand;
                        else
                            tail_nodes[rni].other_side_load = solution.get_route_load_before_included(i) - iDemand;

                        tail_heap.insert(&tail_nodes[rni]);  // heap_insert(rni);
                        rni++;
                        if (rni == max_tail_nodes) {
                            return true;
                        }

                        if constexpr (!direction) {

                            if (iPrev == this->instance.get_depot()) {
                                continue;
                            }

                            tail_nodes[rni] = tail_nodes[rni - 1];
                            tail_nodes[rni].i = iPrev;  // should start from i, but since it is already continued in the previous case, we can skip it
                            tail_nodes[rni].iNext = i;
                            tail_nodes[rni].direction = false;
                            tail_nodes[rni].other_side_load = solution.get_route_load_after_included(j) + iDemand;

                            tail_heap.insert(&tail_nodes[rni]);  // heap_insert(rni);
                            rni++;
                            if (rni == max_tail_nodes) {
                                return true;
                            }
                        }
                    } else if (i_route_feasible) {
                        if (jPrev == this->instance.get_depot()) {
                            continue;
                        }

                        const auto jPrevPrev = solution.get_prev_vertex(jRoute, jPrev);

                        if (jPrevPrev == this->instance.get_depot()) {
                            continue;
                        }
                        tail_nodes[rni].i = jPrev;
                        tail_nodes[rni].iNext = iNext;
                        tail_nodes[rni].direction = false;
                        tail_nodes[rni].infeas_demand = j_route_load - capacity;


                        if constexpr (!direction) tail_nodes[rni].other_side_load = curr.other_side_load - iDemand;
                        else
                            tail_nodes[rni].other_side_load = load_after_iNext;

                        tail_heap.insert(&tail_nodes[rni]);  // heap_insert(rni);
                        ++rni;
                        if (rni == max_tail_nodes) {
                            return true;
                        }

                        if constexpr (direction) {

                            if (iNext == this->instance.get_depot()) {
                                continue;
                            }

                            tail_nodes[rni] = tail_nodes[rni - 1];
                            tail_nodes[rni].i = iNext;
                            tail_nodes[rni].iNext = solution.get_next_vertex(iNext);
                            tail_nodes[rni].direction = true;
                            assert(jPrev != this->instance.get_depot());
                            tail_nodes[rni].other_side_load = solution.get_route_load_before_included(jPrev);

                            tail_heap.insert(&tail_nodes[rni]);  // heap_insert(rni);
                            ++rni;
                            if (rni == max_tail_nodes) {
                                return true;
                            }
                        }
                    }

                    if (delta < -this->tolerance) {
                        curr.other_side_load -= iDemand;
                        return false;
                    }  // new nodes are more promising, get back to the heap and select them
                }

                curr.next_movegen_index = 0;  // now we are changing vertex => consider all its move gen
                if constexpr (direction) {
                    curr.i = curr.iNext;
                } else {
                    curr.i = solution.get_prev_vertex(curr.i);
                }

                if (curr.i == this->instance.get_depot()) {
                    break;
                }

                curr.iNext = solution.get_next_vertex(curr.i);
            }

            // remove curr from heap
            tail_heap.get();  // heap_remove_top();

            return false;
        }


        inline void execute(Solution &solution, __attribute__((unused)) const MoveGenerator &p_move, VertexSet &storage) override {

            // Apply the best chain after storing its vertices
            const auto best_chain_index = feasible_chain;


#ifndef NDEBUG
            auto &best_chain = tail_nodes[best_chain_index];
            auto expected_delta = best_chain.delta_sum;
            auto old_cost = solution.get_cost();
#endif
            // TODO: manca la parte dove si riempie lo storage"

            auto tail_moves = std::vector<int>();
            for (auto ptr = best_chain_index; ptr != -1; ptr = tail_nodes[ptr].predecessor) {
                tail_moves.emplace_back(ptr);
            }

            /* if (tail_moves.size() > 1) {
                store_to_file(this->instance, solution, "./a.sol");
                std::cout << "a\n";
            } */

            for (auto n = static_cast<int>(tail_moves.size()) - 1; n >= 0; --n) {
                const auto ptr = tail_moves[n];
                const auto move = tail_nodes[ptr].move;
                const auto i = move->get_first_vertex();
                const auto j = move->get_second_vertex();

                // std::cout << "Executing tail from MG (" << i << ',' << j << ")\n";

                assert(i != this->instance.get_depot());
                assert(j != this->instance.get_depot());

                const auto iNext = solution.get_next_vertex(i);
                const auto jPrev = solution.get_prev_vertex(j);

                const auto i_route = solution.get_route_index(i);
                const auto j_route = solution.get_route_index(j);

                // std::cout << "Before TAIL\n";
                // solution.print(i_route);
                // solution.print(j_route);

                storage.insert(i);
                storage.insert(iNext);
                storage.insert(jPrev);
                storage.insert(j);
                /* storage.insert(this->instance.get_depot());
                this->update_bits.at(this->instance.get_depot(), UPDATE_BITS_FIRST, true);
                this->update_bits.at(this->instance.get_depot(), UPDATE_BITS_SECOND, true);
                for (auto a = solution.get_first_customer(iRoute); a != this->instance.get_depot(); i = solution.get_next_vertex(a)) {
                    storage.insert(a);
                    this->update_bits.at(a, UPDATE_BITS_FIRST, true);
                    this->update_bits.at(a, UPDATE_BITS_SECOND, true);
                }
                for (auto a = solution.get_first_customer(jRoute); a != this->instance.get_depot(); i = solution.get_next_vertex(a)) {
                    storage.insert(a);
                    this->update_bits.at(a, UPDATE_BITS_FIRST, true);
                    this->update_bits.at(a, UPDATE_BITS_SECOND, true);
                } */

                this->update_bits.at(i, UPDATE_BITS_FIRST, true);
                this->update_bits.at(iNext, UPDATE_BITS_SECOND, true);
                this->update_bits.at(j, UPDATE_BITS_SECOND, true);
                this->update_bits.at(jPrev, UPDATE_BITS_FIRST, true);

                solution.swap_tails(i, i_route, j, j_route);

                // std::cout << "After TAIL\n";
                // solution.print(i_route);
                // solution.print(j_route);
                // std::cout << "\n";

                if (solution.is_route_empty(i_route)) {
                    solution.remove_route(i_route);
                } else {
                    solution.update_cumulative_route_loads(i_route);
                }

                if (solution.is_route_empty(j_route)) {
                    solution.remove_route(j_route);
                } else {
                    solution.update_cumulative_route_loads(j_route);
                }
            }

            assert(std::fabs(old_cost + expected_delta - solution.get_cost()) < 0.01f);

            assert(solution.is_feasible());
            // if (tail_moves.size() > 1) std::cout << "Applied TailChain of size " << tail_moves.size() << "!!!!!\n";
        }

        void post_processing(__attribute__((unused)) Solution &solution) override {
            // REMOVE
            // std::cout << " Called: " << called;
            // std::cout << ", Current FeasChainLenAVG: " << static_cast<double>(feasChainLenSum) / static_cast<double>(foudFeas);
            // std::cout << ", Current feasChainLenMax: " << feasChainLenMax;
            // std::cout << ", Current chainNumSumAVG: " << static_cast<double>(chainNumSum) / static_cast<double>(called);
            // std::cout << std::endl;
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
            int v, next, prev;
            float seq1rem, seq2rem;
        };

        inline Cache12 prepare_cache12(const Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache12();
            c.v = vertex;
            c.next = solution.get_next_vertex(c.v);
            c.prev = solution.get_prev_vertex(c.v);
            c.seq1rem = -solution.get_next_cost(c.v);
            c.seq2rem = -solution.get_prev_cost(c.v);


            return c;
        }

        inline Cache12 prepare_cache12(const Solution &solution, int vertex, int backup) {

            auto c = Cache12();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v, backup);
            c.next = solution.get_next_vertex(route, c.v);
            c.prev = solution.get_prev_vertex(route, c.v);
            c.seq1rem = -solution.get_next_cost(c.v, c.next);
            c.seq2rem = -solution.get_prev_cost(c.v, c.prev);


            return c;
        }

        inline std::pair<float, float> compute_cost_pair(const struct Cache12 i, const struct Cache12 j) {

            const auto c_iv_jv = this->instance.get_cost(i.v, j.v);

            const auto delta1 = i.seq1rem + c_iv_jv + j.seq2rem + this->instance.get_cost(j.prev, i.next);
            const auto delta2 = j.seq1rem + c_iv_jv + i.seq2rem + this->instance.get_cost(i.prev, j.next);

            return {delta1, delta2};
        }

        struct Cache1 {
            int v, next;
            float seq1rem;
        };

        inline Cache1 prepare_cache1(const Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache1();
            c.v = vertex;
            c.next = solution.get_next_vertex(c.v);
            c.seq1rem = -solution.get_next_cost(c.v);

            return c;
        }

        inline Cache1 prepare_cache1(const Solution &solution, int vertex, int backup) {
            auto c = Cache1();
            c.v = vertex;
            const auto route = solution.get_route_index(vertex, backup);
            c.next = solution.get_next_vertex(route, c.v);
            c.seq1rem = -solution.get_next_cost(c.v, c.next);

            return c;
        }

        struct Cache2 {
            int v, prev;
            float seq2rem;
        };

        inline Cache2 prepare_cache2(const Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache2();
            c.v = vertex;
            c.prev = solution.get_prev_vertex(c.v);
            c.seq2rem = -solution.get_prev_cost(c.v);


            return c;
        }

        inline Cache2 prepare_cache2(const Solution &solution, int vertex, int backup) {
            auto c = Cache2();
            c.v = vertex;
            const auto route = solution.get_route_index(vertex, backup);
            c.prev = solution.get_prev_vertex(route, c.v);
            c.seq2rem = -solution.get_prev_cost(c.v, c.prev);


            return c;
        }

        template <typename C1, typename C2>
        inline float compute_cost(const C1 i, const C2 j) {

            const auto c_iv_jv = this->instance.get_cost(i.v, j.v);

            const auto delta1 = i.seq1rem + c_iv_jv + j.seq2rem + this->instance.get_cost(j.prev, i.next);

            return delta1;
        }
    };

}  // namespace cobra

#endif