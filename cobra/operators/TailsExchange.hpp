#ifndef _F4D_TAILSEXCHANGE_HPP_
#define _F4D_TAILSEXCHANGE_HPP_

#include "AbstractOperator.hpp"

namespace cobra {

    template <bool handle_partial_solutions = false>
    class TailsExchange : public AbstractOperator {
    public:
        TailsExchange(const cobra::Instance &instance, MoveGenerators &moves, float tolerance) : AbstractOperator(instance, moves, tolerance) { }

        static constexpr bool is_symmetric = false;

    protected:
        inline void pre_processing(cobra::Solution &solution) override {

            for (int route = solution.get_first_route(); route != cobra::Solution::dummy_route; route = solution.get_next_route(route)) {
                solution.update_cumulative_route_loads(route);
            }
        }

        inline float compute_cost(const cobra::Solution &solution, const MoveGenerator &move) override {
            const auto i = move.get_first_vertex();
            const auto j = move.get_second_vertex();

            const auto iRoute = solution.get_route_index(i, j);
            const auto jRoute = solution.get_route_index(j, i);

            const auto iNext = solution.get_next_vertex(iRoute, i);
            const auto jPrev = solution.get_prev_vertex(jRoute, j);

            const auto delta = -this->instance.get_cost(i, iNext) + this->instance.get_cost(i, j) - this->instance.get_cost(jPrev, j) +
                               this->instance.get_cost(jPrev, iNext) + delta_penalty;

            return delta;
        }

        bool is_feasible(const cobra::Solution &solution, const MoveGenerator &move) override {
            const auto i = move.get_first_vertex();
            const auto j = move.get_second_vertex();

            const auto iRoute = solution.get_route_index(i, j);
            const auto jRoute = solution.get_route_index(j, i);


            return iRoute != jRoute &&
                   // iRoute != jRoute if both i and j are different from the depot
                   solution.get_route_load_before_included(i) + solution.get_route_load_after_included(j) <= this->instance.get_vehicle_capacity() &&
                   solution.get_route_load_before_included(j) - this->instance.get_demand(j) + solution.get_route_load_after_included(i) -
                           this->instance.get_demand(i) <=
                       this->instance.get_vehicle_capacity();
        }

        inline void execute(cobra::Solution &solution, const MoveGenerator &move, cobra::VertexSet &storage) override {

            const auto i = move.get_first_vertex();
            const auto j = move.get_second_vertex();

            const auto iNext = solution.get_next_vertex(i);
            const auto jPrev = solution.get_prev_vertex(j);

            const auto i_route = solution.get_route_index(i);
            const auto j_route = solution.get_route_index(j);

            storage.insert(i);
            storage.insert(iNext);
            storage.insert(jPrev);
            storage.insert(j);

            this->update_bits.at(i, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iNext, UPDATE_BITS_SECOND, true);
            this->update_bits.at(j, UPDATE_BITS_SECOND, true);
            this->update_bits.at(jPrev, UPDATE_BITS_FIRST, true);

            solution.swap_tails(i, i_route, j, j_route);

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

        void post_processing(__attribute__((unused)) cobra::Solution &solution) override { }

        std::string get_additional_statistics() override {
            return std::string();
        }

        struct Cache12 {
            int v, next, prev;
            float seq1rem, seq2rem;
        };

        inline Cache12 prepare_cache12(const cobra::Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache12();
            c.v = vertex;
            c.next = solution.get_next_vertex(c.v);
            c.prev = solution.get_prev_vertex(c.v);
            c.seq1rem = -this->instance.get_cost(c.v, c.next);
            c.seq2rem = -this->instance.get_cost(c.prev, c.v);

            return c;
        }

        inline Cache12 prepare_cache12(const cobra::Solution &solution, int vertex, int backup) {

            auto c = Cache12();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v, backup);
            c.next = solution.get_next_vertex(route, c.v);
            c.prev = solution.get_prev_vertex(route, c.v);
            c.seq1rem = -this->instance.get_cost(c.v, c.next);
            c.seq2rem = -this->instance.get_cost(c.prev, c.v);

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

        inline Cache1 prepare_cache1(const cobra::Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache1();
            c.v = vertex;
            c.next = solution.get_next_vertex(c.v);
            c.seq1rem = -this->instance.get_cost(c.v, c.next);

            return c;
        }

        inline Cache1 prepare_cache1(const cobra::Solution &solution, int vertex, int backup) {
            auto c = Cache1();
            c.v = vertex;
            const auto route = solution.get_route_index(vertex, backup);
            c.next = solution.get_next_vertex(route, c.v);
            c.seq1rem = -this->instance.get_cost(c.v, c.next);

            return c;
        }

        struct Cache2 {
            int v, prev;
            float seq2rem;
        };

        inline Cache2 prepare_cache2(const cobra::Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache2();
            c.v = vertex;
            c.prev = solution.get_prev_vertex(c.v);
            c.seq2rem = -this->instance.get_cost(c.prev, c.v);

            return c;
        }

        inline Cache2 prepare_cache2(const cobra::Solution &solution, int vertex, int backup) {
            auto c = Cache2();
            c.v = vertex;
            const auto route = solution.get_route_index(vertex, backup);
            c.prev = solution.get_prev_vertex(route, c.v);
            c.seq2rem = -this->instance.get_cost(c.prev, c.v);

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