#ifndef _F4D_EJECTIONCHAIN_HPP_
#define _F4D_EJECTIONCHAIN_HPP_

#include "../BinaryHeap.hpp"
#include "../BitMatrix.hpp"
#include "../SmallFlatMap.hpp"
#include "AbstractOperator.hpp"

namespace cobra {

    template <bool log_statistics>
    struct EjectionChainStatistics { };
    template <>
    struct EjectionChainStatistics<true> {

        Welford num_tree_nodes;
        Welford max_tree_depth;
    };

    template <bool handle_partial_solution = false, bool log_statistics = false, int max_relocation_nodes = 100>
    class EjectionChain : EjectionChainStatistics<log_statistics>, public AbstractOperator {

        static const int heap_unheaped = -1;
        static constexpr auto max_chain_length = 5;

       

        struct Relocation {
            short heap_index = heap_unheaped;
            short predecessor = 0;
            int level = 0;
            float delta_sum = 0.0f;
            const MoveGenerator *move = nullptr;
            SmallFlatMap<int, int, 0, 25> modified_routes_loads;
        };

        BitMatrix<2 * max_chain_length + 3> forbidden_i;
        BitMatrix<3 * max_chain_length> forbidden_j;

        std::vector<Relocation> relocation_nodes;
        std::vector<int> feasible_chains;

    
        BinaryHeapPtr<Relocation, &Relocation::heap_index, &Relocation::delta_sum> relo_heap;
        inline int index_of(std::vector<Relocation> &nodes, Relocation &node) {
            return &node - nodes.data();
        }

    public:
        EjectionChain(const Instance &instance_, MoveGenerators &moves_, float tolerance_)
            : AbstractOperator(instance_, moves_, tolerance_),
              forbidden_i(max_relocation_nodes),
              forbidden_j(max_relocation_nodes) {

            relocation_nodes.resize(max_relocation_nodes);

            feasible_chains.reserve(max_relocation_nodes);
           }

        static constexpr bool is_symmetric = false;

    protected:
        inline void pre_processing(__attribute__((unused)) Solution &solution) override { }

        inline float compute_cost(const Solution &solution, const MoveGenerator &move) override {

            const auto i = move.get_first_vertex();
            const auto j = move.get_second_vertex();

            const auto iRoute = solution.get_route_index(i, j);
            const auto jRoute = solution.get_route_index(j, i);

            const auto iPrev = solution.get_prev_vertex(iRoute, i);
            const auto iNext = solution.get_next_vertex(iRoute, i);

            const auto jPrev = solution.get_prev_vertex(jRoute, j);

            auto delta = 0.0f;

            if (j != iNext) {
                delta = -this->instance.get_cost(iPrev, i) - this->instance.get_cost(i, iNext) + this->instance.get_cost(iPrev, iNext) -
                        this->instance.get_cost(jPrev, j) + this->instance.get_cost(jPrev, i) + this->instance.get_cost(i, j);
            }


            return delta;
        }

        // The feasibility step for the ejection chain is more complex than
        // that of other operators. In this context, we use the move generator
        // as the starting point from which a tree of relocations.
        bool is_feasible(const Solution &solution, const MoveGenerator &generating_move) override {

            short rni = 0;  // relocate node index, no. of generated nodes

            feasible_chains.clear();

            // First check whether the current `generating_move` is itself a
            // feasible relocate move. If this is the case, apply it without
            // further searching for feasible ejection chains.

            // The following code is into a block to avoid hiding subsequent
            // variables!
            {

                auto i = generating_move.get_first_vertex();
                auto j = generating_move.get_second_vertex();

                auto iRoute = solution.get_route_index(i, j);
                auto jRoute = solution.get_route_index(j, i);

                auto iPrev = solution.get_prev_vertex(iRoute, i);
                auto iNext = solution.get_next_vertex(iRoute, i);
                auto jPrev = solution.get_prev_vertex(jRoute, j);

                relocation_nodes[rni].move = &generating_move;

                if (iRoute == jRoute || (solution.get_route_load(jRoute) + this->instance.get_demand(i) <= this->instance.get_vehicle_capacity())) {
                    feasible_chains.push_back(0);
                    relocation_nodes[0].predecessor = -1;
                    forbidden_i.reset(0);
                    forbidden_i.set(0, iPrev);
                    forbidden_i.set(0, i);
                    forbidden_i.set(0, iNext);
                    forbidden_i.set(0, jPrev);
                    forbidden_i.set(0, j);
                    return true;
                }

                // If the generating move is not feasible by itself, we start
                // a relocation chain


                // set up state variables
                relocation_nodes[rni].delta_sum = generating_move.get_delta();
                relocation_nodes[rni].level = 0;

                forbidden_i.reset(rni);
                forbidden_i.set(rni, iPrev);
                forbidden_i.set(rni, jPrev);

                forbidden_j.reset(rni);
                forbidden_j.set(rni, i);
                forbidden_j.set(rni, iNext);
                forbidden_j.set(rni, j);

                relocation_nodes[rni].modified_routes_loads.clear();
                relocation_nodes[rni].modified_routes_loads[iRoute] = solution.get_route_load(iRoute) - this->instance.get_demand(i);
                relocation_nodes[rni].modified_routes_loads[jRoute] = solution.get_route_load(jRoute) + this->instance.get_demand(i);
                relocation_nodes[rni].predecessor = -1;

                relo_heap.reset();                         
                relo_heap.insert(&relocation_nodes[rni]); 
                rni++;
            }

            while (!relo_heap.empty()) {

                auto &curr = *relo_heap.get();                             
                const auto curr_index = index_of(relocation_nodes, curr);  

                // retrieve the route from which we would like to remove some vertex
                const auto iRoute = solution.get_route_index(curr.move->get_second_vertex());

                // retrieve the updated 'iRoute' load (iRoute will always be in the map)
                assert(curr.modified_routes_loads.count(iRoute));
                const auto iRoute_load = curr.modified_routes_loads[iRoute];

                // scan 'iRoute' searching for customers to remove, which will make the route feasible
                for (auto i = solution.get_first_customer(iRoute); i != this->instance.get_depot(); i = solution.get_next_vertex(i)) {

                    // check whether removing 'i' is sufficient to restore the 'iRoute' feasibility
                    const auto iDemand = this->instance.get_demand(i);
                    if (iRoute_load - iDemand > this->instance.get_vehicle_capacity()) {
                        continue;
                    }

                    // check whether this vertex can be used as 'i'
                    if (forbidden_i.is_set(curr_index, i) || forbidden_j.is_set(curr_index, i)) {
                        continue;
                    }

                    // compute once:
                    const auto iPrev = solution.get_prev_vertex(iRoute, i);
                    const auto iNext = solution.get_next_vertex(iRoute, i);
                    const auto iCost = -this->instance.get_cost(iPrev, i) - this->instance.get_cost(i, iNext) + this->instance.get_cost(iPrev, iNext);
                    // retrieve the available generating_move generators involving 'i'. They may be {i j} and {j i}
                    // scan them ...
                    for (const auto move_index : this->moves.get_move_generator_indices_involving_1st(i)) {
                        auto &move = this->moves.get(move_index);

                        if constexpr (handle_partial_solution) {
                            if (!solution.is_vertex_in_solution(move.get_first_vertex()) || !solution.is_vertex_in_solution(move.get_second_vertex())) {
                                continue;
                            }
                        }

                        assert(move.get_first_vertex() == i);


                        const auto j = move.get_second_vertex();

                        // check whether this vertex can be used as 'j', we cannot use the depot because we want to relocate in another route
                        if (j == this->instance.get_depot() || forbidden_j.is_set(curr_index, j)) {
                            continue;
                        }

                        const auto jRoute = solution.get_route_index(j);

                        // relocate in another route! Note we are removing 'i' to restore 'iRoute' feasibility
                        if (jRoute == iRoute) {
                            continue;
                        }

                        // ========================================================================
                        // For move generators not in the 'MovesHeap' we don't know whether
                        // their delta was updated during the initialization stage or not. To this end,
                        // we just recompute the value of 'delta'. Note that previous checks ensure
                        // to be in a particular scenario, i.e. we know iRoute != jRoute and this
                        // simplifies the computation of 'delta'.
                        //
                        // The ejection chain operation is the only move affected by this problem
                        // because it is a combination of moves all but the generating one that
                        // might not be in the 'MovesHeap'.
                        // ========================================================================
                        const auto jPrev = solution.get_prev_vertex(jRoute, j);
                        const auto correct_delta = iCost - this->instance.get_cost(jPrev, j) + this->instance.get_cost(jPrev, i) +
                                                   this->instance.get_cost(i, j);
                        move.set_delta(correct_delta);

                        // ... that keep the chain improving
                        if (move.get_delta() + curr.delta_sum > -this->tolerance) {
                            continue;
                        }

                        // check whether we have already worked with 'jRoute'
                        auto jRoute_load = [&] {
                            if (const auto &pair = curr.modified_routes_loads.find(jRoute); pair.first != 0) {
                                return pair.second;
                            } else {
                                return solution.get_route_load(jRoute);
                            }
                        }();

                        assert(std::abs(move.get_delta() - this->compute_cost(solution, move)) < 0.01f);

                        bool feas = jRoute_load + iDemand <= this->instance.get_vehicle_capacity();

                        if (feas || curr.level + 1 < max_chain_length - 1) {
                            // at this point, the move {i j} is able to restore 'iRoute' feasibility, store it in 'rni' position
                            relocation_nodes[rni].move = &move;
                            relocation_nodes[rni].delta_sum = curr.delta_sum + move.get_delta();
                            relocation_nodes[rni].level = curr.level + 1;

                            forbidden_i.overwrite(curr_index, rni);
                            forbidden_i.set(rni, iPrev);
                            forbidden_i.set(rni, jPrev);

                            forbidden_j.overwrite(curr_index, rni);
                            forbidden_j.set(rni, i);
                            forbidden_j.set(rni, iNext);
                            forbidden_j.set(rni, j);

                            relocation_nodes[rni].modified_routes_loads = curr.modified_routes_loads;
                            relocation_nodes[rni].modified_routes_loads[iRoute] = iRoute_load - iDemand;
                            relocation_nodes[rni].modified_routes_loads[jRoute] = jRoute_load + iDemand;

                            relocation_nodes[rni].predecessor = curr_index;
                            relo_heap.insert(&relocation_nodes[rni]);  

                            // if also 'jRoute' is feasible we have found a feasible chain!
                            if (feas) {
                                feasible_chains.push_back(rni);
                                goto end;
                            }
                        }

                        rni++;

                        if (rni == max_relocation_nodes) {
                            goto end;
                        }
                    }
                }
            }

        end:

            if constexpr (log_statistics) {

                auto max_depth = 0;
                for (auto idx = 0; idx < rni; idx++) {
                    if (relocation_nodes[idx].level > max_depth) {
                        max_depth = relocation_nodes[idx].level;
                    }
                }
                this->max_tree_depth.update(max_depth);

                this->num_tree_nodes.update(rni);
            }


            return !feasible_chains.empty();
        }

        inline void execute(Solution &solution, __attribute__((unused)) const MoveGenerator &p_move, VertexSet &storage) override {


            // Apply the best chain after storing its vertices
            const auto best_chain_index = feasible_chains[0];

            for (auto i : forbidden_i.get_set_entries_possibly_with_duplicates(best_chain_index)) {
                storage.insert(i);
            }
            for (auto j : forbidden_j.get_set_entries_possibly_with_duplicates(best_chain_index)) {
                storage.insert(j);
            }

            for (auto ptr = best_chain_index; ptr != -1; ptr = relocation_nodes[ptr].predecessor) {
                const auto &move = relocation_nodes[ptr].move;

                const auto i = move->get_first_vertex();
                const auto j = move->get_second_vertex();

                const auto iRoute = solution.get_route_index(i, j);
                const auto jRoute = solution.get_route_index(j, i);


                this->update_bits.at(solution.get_prev_vertex(iRoute, i), UPDATE_BITS_FIRST, true);
                this->update_bits.at(i, UPDATE_BITS_FIRST, true);
                this->update_bits.at(i, UPDATE_BITS_SECOND, true);
                const auto iNext = solution.get_next_vertex(iRoute, i);
                this->update_bits.at(iNext, UPDATE_BITS_FIRST, true);
                this->update_bits.at(iNext, UPDATE_BITS_SECOND, true);
                this->update_bits.at(j, UPDATE_BITS_FIRST, true);
                this->update_bits.at(j, UPDATE_BITS_SECOND, true);
                this->update_bits.at(solution.get_prev_vertex(jRoute, j), UPDATE_BITS_FIRST, true);

                solution.remove_vertex(iRoute, i);
                solution.insert_vertex_before(jRoute, j, i);

                if (solution.is_route_empty(iRoute)) {
                    solution.remove_route(iRoute);
                }
            }

            assert(solution.is_feasible());
        }

        void post_processing(__attribute__((unused)) Solution &solution) override {
            
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
            int v, prev, next;
            float vrem, prevrem;
        };

        inline Cache12 prepare_cache12(const Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache12();
            c.v = vertex;
            c.prev = solution.get_prev_vertex(c.v);
            c.next = solution.get_next_vertex(c.v);

            // c.v = i in (i, j)
            c.vrem = -this->instance.get_cost(c.prev, c.v) - this->instance.get_cost(c.v, c.next) + this->instance.get_cost(c.prev, c.next);
            // c.v = j in (i, j)
            c.prevrem = -this->instance.get_cost(c.prev, c.v);

            return c;
        }


        inline Cache12 prepare_cache12(const Solution &solution, int vertex, int backup) {

            auto c = Cache12();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v, backup);
            c.prev = solution.get_prev_vertex(route, c.v);
            c.next = solution.get_next_vertex(route, c.v);

            // c.v = i in (i, j)
            c.vrem = -this->instance.get_cost(c.prev, c.v) - this->instance.get_cost(c.v, c.next) + this->instance.get_cost(c.prev, c.next);
            // c.v = j in (i, j)
            c.prevrem = -this->instance.get_cost(c.prev, c.v);

            return c;
        }

        inline std::pair<float, float> compute_cost_pair(const struct Cache12 i, const struct Cache12 j) {

            const auto c_iv_jv = this->instance.get_cost(i.v, j.v);
            const auto delta1 = j.v != i.next ? i.vrem + j.prevrem + this->instance.get_cost(j.prev, i.v) + c_iv_jv : 0.0f;
            const auto delta2 = i.v != j.next ? j.vrem + i.prevrem + this->instance.get_cost(i.prev, j.v) + c_iv_jv : 0.0f;

            return {delta1, delta2};
        }

        struct Cache1 {
            int v, prev, next;
            float vrem;
        };

        inline Cache1 prepare_cache1(const Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache1();
            c.v = vertex;
            c.prev = solution.get_prev_vertex(c.v);
            c.next = solution.get_next_vertex(c.v);
            c.vrem = -this->instance.get_cost(c.prev, c.v) - this->instance.get_cost(c.v, c.next) + this->instance.get_cost(c.prev, c.next);
            return c;
        }

        inline Cache1 prepare_cache1(const Solution &solution, int vertex, int backup) {
            auto c = Cache1();
            c.v = vertex;
            const auto route = solution.get_route_index(vertex, backup);
            c.prev = solution.get_prev_vertex(route, c.v);
            c.next = solution.get_next_vertex(route, c.v);
            c.vrem = -this->instance.get_cost(c.prev, c.v) - this->instance.get_cost(c.v, c.next) + this->instance.get_cost(c.prev, c.next);
            return c;
        }

        struct Cache2 {
            int v, prev;
            float prevrem;
        };

        inline Cache2 prepare_cache2(const Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache2();
            c.v = vertex;
            c.prev = solution.get_prev_vertex(c.v);
            c.prevrem = -this->instance.get_cost(c.prev, c.v);
            return c;
        }

        inline Cache2 prepare_cache2(const Solution &solution, int vertex, int backup) {
            auto c = Cache2();
            c.v = vertex;
            const auto route = solution.get_route_index(vertex, backup);
            c.prev = solution.get_prev_vertex(route, c.v);
            c.prevrem = -this->instance.get_cost(c.prev, c.v);
            return c;
        }

        template <typename C1, typename C2>
        inline float compute_cost(const C1 i, const C2 j) {

            return j.v != i.next ? i.vrem + j.prevrem + this->instance.get_cost(j.prev, i.v) + this->instance.get_cost(i.v, j.v) : 0.0f;
        }
    };

}  // namespace cobra

#endif