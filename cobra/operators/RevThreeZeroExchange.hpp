#ifndef _F4D_REVTHREEZEROEXCHANGE_HPP_
#define _F4D_REVTHREEZEROEXCHANGE_HPP_

#include "AbstractOperator.hpp"

namespace cobra {

    template <bool handle_partial_solutions = false>
    class RevThreeZeroExchange : public AbstractOperator {
    public:
        RevThreeZeroExchange(const cobra::Instance &instance, MoveGenerators &moves, float tolerance) : AbstractOperator(instance, moves, tolerance) { }

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
            const auto iPrevPrevPrev = solution.get_prev_vertex(iRoute, iPrevPrev);

            const auto jNext = solution.get_next_vertex(jRoute, j);

            const auto iSequenceRem = -this->instance.get_cost(iPrevPrevPrev, iPrevPrev) - this->instance.get_cost(i, iNext);
            const auto jSequenceRem = -this->instance.get_cost(j, jNext);
            const auto iSequenceAdd = +this->instance.get_cost(jNext, iPrevPrev) + this->instance.get_cost(i, j);
            const auto iFilling = +this->instance.get_cost(iPrevPrevPrev, iNext);

            return iSequenceAdd + iFilling + iSequenceRem + jSequenceRem;
        }

        bool is_feasible(const cobra::Solution &solution, const MoveGenerator &move) override {
            const auto i = move.get_first_vertex();
            const auto j = move.get_second_vertex();

            const auto iRoute = solution.get_route_index(i, j);
            const auto jRoute = solution.get_route_index(j, i);

            const auto iPrev = solution.get_prev_vertex(iRoute, i);
            const auto iPrevPrev = solution.get_prev_vertex(iRoute, iPrev);


            return (iRoute != jRoute && iPrev != this->instance.get_depot() && iPrevPrev != this->instance.get_depot() &&
                    solution.get_route_load(jRoute) + this->instance.get_demand(i) + this->instance.get_demand(iPrev) + this->instance.get_demand(iPrevPrev) <=
                        this->instance.get_vehicle_capacity()) ||
                   (iRoute == jRoute && j != iPrev && j != iPrevPrev && j != solution.get_prev_vertex(iRoute, iPrevPrev));
        }

        inline void execute(cobra::Solution &solution, const MoveGenerator &move, cobra::VertexSet &storage) override {
            const auto i = move.get_first_vertex();
            const auto j = move.get_second_vertex();

            const auto iRoute = solution.get_route_index(i, j);
            const auto jRoute = solution.get_route_index(j, i);

            const auto iPrev = solution.get_prev_vertex(iRoute, i);
            const auto iPrevPrev = solution.get_prev_vertex(iRoute, iPrev);
            const auto iPrevPrevPrev = solution.get_prev_vertex(iRoute, iPrevPrev);

            const auto iNext = solution.get_next_vertex(iRoute, i);
            const auto iNextNext = solution.get_next_vertex(iRoute, iNext);
            const auto iNextNextNext = solution.get_next_vertex(iRoute, iNextNext);

            const auto jNext = solution.get_next_vertex(jRoute, j);
            const auto jNextNext = solution.get_next_vertex(jRoute, jNext);
            const auto jNextNextNext = solution.get_next_vertex(jRoute, jNextNext);


            storage.insert(iPrevPrevPrev);
            storage.insert(iPrevPrev);
            storage.insert(iPrev);
            storage.insert(i);
            storage.insert(iNext);
            storage.insert(iNextNext);
            storage.insert(iNextNextNext);
            storage.insert(j);
            storage.insert(jNext);
            storage.insert(jNextNext);
            storage.insert(jNextNextNext);

            this->update_bits.at(iPrevPrevPrev, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iPrevPrevPrev, UPDATE_BITS_SECOND, true);
            this->update_bits.at(iPrevPrev, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iPrevPrev, UPDATE_BITS_SECOND, true);  // because of reversal
            this->update_bits.at(iPrev, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iPrev, UPDATE_BITS_SECOND, true);  // because of reversal
            this->update_bits.at(i, UPDATE_BITS_FIRST, true);
            this->update_bits.at(i, UPDATE_BITS_SECOND, true);
            this->update_bits.at(iNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iNextNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iNextNextNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(jNextNextNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(jNextNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(jNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(j, UPDATE_BITS_FIRST, true);
            this->update_bits.at(j, UPDATE_BITS_SECOND, true);


            solution.remove_vertex(iRoute, i);
            solution.remove_vertex(iRoute, iPrev);
            solution.remove_vertex(iRoute, iPrevPrev);

            solution.insert_vertex_before(jRoute, jNext, i);
            solution.insert_vertex_before(jRoute, jNext, iPrev);
            solution.insert_vertex_before(jRoute, jNext, iPrevPrev);

            if (solution.is_route_empty(iRoute)) {
                solution.remove_route(iRoute);
            }
        }

        void post_processing(__attribute__((unused)) cobra::Solution &solution) override { }

        std::string get_additional_statistics() override {
            return std::string();
        }

        struct Cache12 {
            int v, next, prevprev, prevprevprev;
            float sequrem, nextrem;
        };

        inline Cache12 prepare_cache12(const cobra::Solution &solution, int vertex) {

            assert(vertex != this->instance.get_depot());
            auto c = Cache12();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v);
            const auto prev = solution.get_prev_vertex(c.v);
            c.next = solution.get_next_vertex(c.v);
            c.prevprev = solution.get_prev_vertex(route, prev);
            c.prevprevprev = solution.get_prev_vertex(route, c.prevprev);

            const auto c_v_next = this->instance.get_cost(c.v, c.next);

            // c.v = i in (i, j)
            c.sequrem = -this->instance.get_cost(c.prevprevprev, c.prevprev) - c_v_next + this->instance.get_cost(c.prevprevprev, c.next);
            // c.v = j in (i, j)
            c.nextrem = -c_v_next;

            return c;
        }


        inline Cache12 prepare_cache12(const cobra::Solution &solution, int vertex, int backup) {

            auto c = Cache12();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v, backup);
            const auto prev = solution.get_prev_vertex(route, c.v);
            c.next = solution.get_next_vertex(route, c.v);
            c.prevprev = solution.get_prev_vertex(route, prev);
            c.prevprevprev = solution.get_prev_vertex(route, c.prevprev);

            const auto c_v_next = this->instance.get_cost(c.v, c.next);

            // c.v = i in (i, j)
            c.sequrem = -this->instance.get_cost(c.prevprevprev, c.prevprev) - c_v_next + this->instance.get_cost(c.prevprevprev, c.next);
            // c.v = j in (i, j)
            c.nextrem = -c_v_next;

            return c;
        }

        inline std::pair<float, float> compute_cost_pair(const struct Cache12 i, const struct Cache12 j) {

            const auto c_iv_jv = this->instance.get_cost(i.v, j.v);

            const auto delta1 = this->instance.get_cost(j.next, i.prevprev) + c_iv_jv + i.sequrem + j.nextrem;
            const auto delta2 = this->instance.get_cost(i.next, j.prevprev) + c_iv_jv + j.sequrem + i.nextrem;

            return {delta1, delta2};
        }

        struct Cache1 {
            int v, prevprev;
            float sequrem;
        };

        inline Cache1 prepare_cache1(const cobra::Solution &solution, int vertex) {

            assert(vertex != this->instance.get_depot());
            auto c = Cache1();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v);
            const auto prev = solution.get_prev_vertex(c.v);
            const auto next = solution.get_next_vertex(c.v);
            c.prevprev = solution.get_prev_vertex(route, prev);
            const auto prevprevprev = solution.get_prev_vertex(route, c.prevprev);

            c.sequrem = -this->instance.get_cost(prevprevprev, c.prevprev) - this->instance.get_cost(c.v, next) + this->instance.get_cost(prevprevprev, next);

            return c;
        }

        inline Cache1 prepare_cache1(const cobra::Solution &solution, int vertex, int backup) {

            auto c = Cache1();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v, backup);
            const auto prev = solution.get_prev_vertex(route, c.v);
            const auto next = solution.get_next_vertex(route, c.v);
            c.prevprev = solution.get_prev_vertex(route, prev);
            const auto prevprevprev = solution.get_prev_vertex(route, c.prevprev);

            c.sequrem = -this->instance.get_cost(prevprevprev, c.prevprev) - this->instance.get_cost(c.v, next) + this->instance.get_cost(prevprevprev, next);

            return c;
        }

        struct Cache2 {
            int v, next;
            float nextrem;
        };

        inline Cache2 prepare_cache2(const cobra::Solution &solution, int vertex) {

            assert(vertex != this->instance.get_depot());
            auto c = Cache2();
            c.v = vertex;
            c.next = solution.get_next_vertex(c.v);

            c.nextrem = -this->instance.get_cost(c.v, c.next);

            return c;
        }

        inline Cache2 prepare_cache2(const cobra::Solution &solution, int vertex, int backup) {

            auto c = Cache2();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v, backup);
            c.next = solution.get_next_vertex(route, c.v);

            c.nextrem = -this->instance.get_cost(c.v, c.next);

            return c;
        }

        template <typename C1, typename C2>
        inline float compute_cost(const C1 i, const C2 j) {

            const auto c_iv_jv = this->instance.get_cost(i.v, j.v);
            const auto delta = this->instance.get_cost(j.next, i.prevprev) + c_iv_jv + i.sequrem + j.nextrem;

            return delta;
        }
    };

}  // namespace cobra

#endif