#ifndef _F4D_REVTWOZEROEXCHANGE_HPP_
#define _F4D_REVTWOZEROEXCHANGE_HPP_

#include "AbstractOperator.hpp"

namespace cobra {

    template <bool handle_partial_solutions = false>
    class RevTwoZeroExchange : public AbstractOperator {
    public:
        RevTwoZeroExchange(const cobra::Instance &instance, MoveGenerators &moves, float tolerance) : AbstractOperator(instance, moves, tolerance) { }

        static constexpr bool is_symmetric = false;

    protected:
        inline void pre_processing(__attribute__((unused)) cobra::Solution &solution) override { }

        inline float compute_cost(const cobra::Solution &solution, const MoveGenerator &move) override {
            const auto i = move.get_first_vertex();
            const auto j = move.get_second_vertex();

            const auto iRoute = solution.get_route_index(i, j);
            const auto jRoute = solution.get_route_index(j, i);

            const auto iNext = solution.get_next_vertex(iRoute, i);
            const auto iPrev = solution.get_prev_vertex(iRoute, i);
            const auto iPrevPrev = solution.get_prev_vertex(iRoute, iPrev);

            const auto jNext = solution.get_next_vertex(jRoute, j);

            return -this->instance.get_cost(iPrevPrev, iPrev) - this->instance.get_cost(i, iNext) + this->instance.get_cost(iPrevPrev, iNext) -
                   this->instance.get_cost(j, jNext) + this->instance.get_cost(i, j) + this->instance.get_cost(iPrev, jNext);
        }

        bool is_feasible(const cobra::Solution &solution, const MoveGenerator &move) override {
            const auto i = move.get_first_vertex();
            const auto j = move.get_second_vertex();

            const auto iRoute = solution.get_route_index(i, j);
            const auto jRoute = solution.get_route_index(j, i);

            const auto iPrev = solution.get_prev_vertex(iRoute, i);

            return (iRoute != jRoute && iPrev != this->instance.get_depot() &&
                    solution.get_route_load(jRoute) + this->instance.get_demand(i) + this->instance.get_demand(iPrev) <=
                        this->instance.get_vehicle_capacity()) ||
                   (iRoute == jRoute && iPrev != j && j != solution.get_prev_vertex(iRoute, iPrev));
        }

        inline void execute(cobra::Solution &solution, const MoveGenerator &move, cobra::VertexSet &storage) override {

            const auto i = move.get_first_vertex();
            const auto j = move.get_second_vertex();


            const auto iRoute = solution.get_route_index(i, j);

            const auto iPrev = solution.get_prev_vertex(iRoute, i);
            const auto iPrevPrev = solution.get_prev_vertex(iRoute, iPrev);
            const auto iNext = solution.get_next_vertex(iRoute, i);
            const auto iNextNext = solution.get_next_vertex(iRoute, iNext);

            const auto jRoute = solution.get_route_index(j, i);

            const auto jNext = solution.get_next_vertex(jRoute, j);
            const auto jNextNext = solution.get_next_vertex(jRoute, jNext);

            storage.insert(iPrevPrev);
            storage.insert(iPrev);
            storage.insert(i);
            storage.insert(iNext);
            storage.insert(iNextNext);
            storage.insert(j);
            storage.insert(jNext);
            storage.insert(jNextNext);

            this->update_bits.at(iPrevPrev, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iPrevPrev, UPDATE_BITS_SECOND, true);
            this->update_bits.at(iPrev, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iPrev, UPDATE_BITS_SECOND, true);  // predecessor of iPrev changes due to the reversal
            this->update_bits.at(i, UPDATE_BITS_FIRST, true);
            this->update_bits.at(i, UPDATE_BITS_SECOND, true);
            this->update_bits.at(iNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iNextNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(j, UPDATE_BITS_FIRST, true);
            this->update_bits.at(j, UPDATE_BITS_SECOND, true);
            this->update_bits.at(jNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(jNextNext, UPDATE_BITS_FIRST, true);


            solution.remove_vertex(iRoute, iPrev);
            solution.remove_vertex(iRoute, i);
            solution.insert_vertex_before(jRoute, jNext, i);
            solution.insert_vertex_before(jRoute, jNext, iPrev);

            if (solution.is_route_empty(iRoute)) {
                solution.remove_route(iRoute);
            }
        }

        void post_processing(__attribute__((unused)) cobra::Solution &solution) override { }

        std::string get_additional_statistics() override {
            return std::string();
        }


        struct Cache12 {
            int v, prev, next;
            float seqrem, prevrem;
        };

        inline Cache12 prepare_cache12(const cobra::Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache12();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v);
            c.prev = solution.get_prev_vertex(c.v);
            const auto prevprev = solution.get_prev_vertex(route, c.prev);
            c.next = solution.get_next_vertex(c.v);

            c.prevrem = -this->instance.get_cost(c.v, c.next);
            c.seqrem = -this->instance.get_cost(prevprev, c.prev) + c.prevrem + this->instance.get_cost(prevprev, c.next);


            return c;
        }


        inline Cache12 prepare_cache12(const cobra::Solution &solution, int vertex, int backup) {

            auto c = Cache12();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v, backup);
            c.prev = solution.get_prev_vertex(route, c.v);
            const auto prevprev = solution.get_prev_vertex(route, c.prev);
            c.next = solution.get_next_vertex(route, c.v);

            c.prevrem = -this->instance.get_cost(c.v, c.next);
            c.seqrem = -this->instance.get_cost(prevprev, c.prev) + c.prevrem + this->instance.get_cost(prevprev, c.next);


            return c;
        }

        inline std::pair<float, float> compute_cost_pair(const struct Cache12 i, const struct Cache12 j) {

            const auto c_iv_jv = this->instance.get_cost(i.v, j.v);

            const auto iSequenceAdd = c_iv_jv + this->instance.get_cost(i.prev, j.next);
            const auto jSequenceAdd = c_iv_jv + this->instance.get_cost(j.prev, i.next);

            const auto delta1 = iSequenceAdd + i.seqrem + j.prevrem;
            const auto delta2 = jSequenceAdd + j.seqrem + i.prevrem;

            return {delta1, delta2};
        }

        struct Cache1 {
            int v, prev;
            float seqrem;
        };

        inline Cache1 prepare_cache1(const cobra::Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache1();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v);
            c.prev = solution.get_prev_vertex(c.v);
            const auto prevprev = solution.get_prev_vertex(route, c.prev);
            const auto next = solution.get_next_vertex(c.v);

            c.seqrem = -this->instance.get_cost(prevprev, c.prev) - this->instance.get_cost(c.v, next) + this->instance.get_cost(prevprev, next);


            return c;
        }

        inline Cache1 prepare_cache1(const cobra::Solution &solution, int vertex, int backup) {
            auto c = Cache1();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v, backup);
            c.prev = solution.get_prev_vertex(route, c.v);
            const auto prevprev = solution.get_prev_vertex(route, c.prev);
            const auto next = solution.get_next_vertex(route, c.v);

            c.seqrem = -this->instance.get_cost(prevprev, c.prev) - this->instance.get_cost(c.v, next) + this->instance.get_cost(prevprev, next);


            return c;
        }

        struct Cache2 {
            int v, next;
            float prevrem;
        };

        inline Cache2 prepare_cache2(const cobra::Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache2();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v);
            c.next = solution.get_next_vertex(c.v);

            c.prevrem = -this->instance.get_cost(c.v, c.next);


            return c;
        }

        inline Cache2 prepare_cache2(const cobra::Solution &solution, int vertex, int backup) {

            auto c = Cache2();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v, backup);
            c.next = solution.get_next_vertex(route, c.v);

            c.prevrem = -this->instance.get_cost(c.v, c.next);


            return c;
        }


        template <typename C1, typename C2>
        inline float compute_cost(const C1 i, const C2 j) {

            const auto iSequenceAdd = +this->instance.get_cost(i.v, j.v) + this->instance.get_cost(i.prev, j.next);

            const auto delta = iSequenceAdd + i.seqrem + j.prevrem;

            return delta;
        }
    };

}  // namespace cobra

#endif