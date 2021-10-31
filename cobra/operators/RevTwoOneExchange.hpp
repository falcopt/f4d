#ifndef _F4D_REVTWOONEEXCHANGE_HPP_
#define _F4D_REVTWOONEEXCHANGE_HPP_

#include "AbstractOperator.hpp"

namespace cobra {

    template <bool handle_partial_solutions = false>
    class RevTwoOneExchange : public AbstractOperator {
    public:
        RevTwoOneExchange(const cobra::Instance &instance, MoveGenerators &moves, float tolerance) : AbstractOperator(instance, moves, tolerance) { }

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
            const auto jNextNext = solution.get_next_vertex(jRoute, jNext);

            const auto iSequenceRem = -this->instance.get_cost(iPrevPrev, iPrev) - this->instance.get_cost(i, iNext);
            const auto jPrevRem = -this->instance.get_cost(j, jNext) - this->instance.get_cost(jNext, jNextNext);

            const auto iSequenceAdd = +this->instance.get_cost(jNextNext, iPrev) + this->instance.get_cost(i, j);
            const auto jNextAdd = +this->instance.get_cost(iPrevPrev, jNext) + this->instance.get_cost(jNext, iNext);

            const auto delta = iSequenceAdd + jNextAdd + iSequenceRem + jPrevRem + delta_penalty;

            return delta;
        }

        bool is_feasible(const cobra::Solution &solution, const MoveGenerator &move) override {
            const auto i = move.get_first_vertex();
            const auto j = move.get_second_vertex();

            const auto iRoute = solution.get_route_index(i, j);
            const auto jRoute = solution.get_route_index(j, i);

            const auto iPrev = solution.get_prev_vertex(iRoute, i);
            const auto iPrevPrev = solution.get_prev_vertex(iRoute, iPrev);
            const auto jNext = solution.get_next_vertex(jRoute, j);


            return (iRoute != jRoute && iPrev != this->instance.get_depot() && jNext != this->instance.get_depot() &&
                    solution.get_route_load(jRoute) - this->instance.get_demand(jNext) + this->instance.get_demand(iPrev) + this->instance.get_demand(i) <=
                        this->instance.get_vehicle_capacity() &&
                    solution.get_route_load(iRoute) + this->instance.get_demand(jNext) - this->instance.get_demand(iPrev) - this->instance.get_demand(i) <=
                        this->instance.get_vehicle_capacity()) ||
                   (iRoute == jRoute && j != iPrev && j != iPrevPrev && jNext != iPrevPrev);
        }

        inline void execute(cobra::Solution &solution, const MoveGenerator &move, cobra::VertexSet &storage) override {

            const auto i = move.get_first_vertex();
            const auto j = move.get_second_vertex();

            const auto iRoute = solution.get_route_index(i, j);
            const auto jRoute = solution.get_route_index(j, i);

            const auto iNext = solution.get_next_vertex(iRoute, i);
            const auto iPrev = solution.get_prev_vertex(iRoute, i);
            const auto iPrevPrev = solution.get_prev_vertex(iRoute, iPrev);
            const auto iPrevPrevPrev = solution.get_prev_vertex(iRoute, iPrevPrev);

            const auto iNextNext = solution.get_next_vertex(iRoute, iNext);

            const auto jPrev = solution.get_prev_vertex(jRoute, j);
            const auto jNext = solution.get_next_vertex(jRoute, j);
            const auto jNextNext = solution.get_next_vertex(jRoute, jNext);
            const auto jNextNextNext = solution.get_next_vertex(jRoute, jNextNext);


            storage.insert(iPrevPrevPrev);
            storage.insert(iPrevPrev);
            storage.insert(iPrev);
            storage.insert(i);
            storage.insert(iNext);
            storage.insert(iNextNext);
            storage.insert(jPrev);
            storage.insert(j);
            storage.insert(jNext);
            storage.insert(jNextNext);
            storage.insert(jNextNextNext);

            this->update_bits.at(iPrevPrevPrev, UPDATE_BITS_SECOND, true);
            this->update_bits.at(iPrevPrev, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iPrevPrev, UPDATE_BITS_SECOND, true);
            this->update_bits.at(iPrev, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iPrev, UPDATE_BITS_SECOND, true);
            this->update_bits.at(i, UPDATE_BITS_FIRST, true);
            this->update_bits.at(i, UPDATE_BITS_SECOND, true);
            this->update_bits.at(iNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iNextNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(jNextNextNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(jNextNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(jNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(jNext, UPDATE_BITS_SECOND, true);
            this->update_bits.at(j, UPDATE_BITS_FIRST, true);
            this->update_bits.at(j, UPDATE_BITS_SECOND, true);
            this->update_bits.at(jPrev, UPDATE_BITS_SECOND, true);


            solution.remove_vertex(jRoute, jNext);
            solution.insert_vertex_before(iRoute, iNext, jNext);

            solution.remove_vertex(iRoute, i);
            solution.remove_vertex(iRoute, iPrev);

            solution.insert_vertex_before(jRoute, jNextNext, i);
            solution.insert_vertex_before(jRoute, jNextNext, iPrev);
        }

        void post_processing(__attribute__((unused)) cobra::Solution &solution) override { }

        std::string get_additional_statistics() override {
            return std::string();
        }

        /*** NEW SECTION ***/

        struct Cache12 {
            int v, prev, prevprev, next, nextnext;
            float seqrem, prevrem;
        };

        inline Cache12 prepare_cache12(const cobra::Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache12();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v);
            c.prev = solution.get_prev_vertex(c.v);
            c.prevprev = solution.get_prev_vertex(route, c.prev);
            c.next = solution.get_next_vertex(c.v);
            c.nextnext = solution.get_next_vertex(route, c.next);

            const auto c_v_next = this->instance.get_cost(c.v, c.next);
            c.seqrem = -this->instance.get_cost(c.prevprev, c.prev) - c_v_next;
            c.prevrem = -c_v_next - this->instance.get_cost(c.next, c.nextnext);

            return c;
        }

        inline Cache12 prepare_cache12(const cobra::Solution &solution, int vertex, int backup) {

            auto c = Cache12();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v, backup);
            c.prev = solution.get_prev_vertex(route, c.v);
            c.prevprev = solution.get_prev_vertex(route, c.prev);
            c.next = solution.get_next_vertex(route, c.v);
            c.nextnext = solution.get_next_vertex(route, c.next);

            const auto c_v_next = this->instance.get_cost(c.v, c.next);
            c.seqrem = -this->instance.get_cost(c.prevprev, c.prev) - c_v_next;
            c.prevrem = -c_v_next - this->instance.get_cost(c.next, c.nextnext);

            return c;
        }

        inline std::pair<float, float> compute_cost_pair(const struct Cache12 i, const struct Cache12 j) {

            const auto c_iv_jv = this->instance.get_cost(i.v, j.v);
            const auto c_inext_jnext = this->instance.get_cost(i.next, j.next);

            const auto iSequenceAdd = +this->instance.get_cost(j.nextnext, i.prev) + c_iv_jv;
            const auto jSequenceAdd = +this->instance.get_cost(i.nextnext, j.prev) + c_iv_jv;

            const auto jNextAdd = +this->instance.get_cost(i.prevprev, j.next) + c_inext_jnext;
            const auto iNextAdd = +this->instance.get_cost(j.prevprev, i.next) + c_inext_jnext;


            const auto delta1 = iSequenceAdd + jNextAdd + i.seqrem + j.prevrem;
            const auto delta2 = jSequenceAdd + iNextAdd + j.seqrem + i.prevrem;

            return {delta1, delta2};
        }

        struct Cache1 {
            int v, prev, prevprev, next;
            float seqrem;
        };

        inline Cache1 prepare_cache1(const cobra::Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache1();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v);
            c.prev = solution.get_prev_vertex(c.v);
            c.prevprev = solution.get_prev_vertex(route, c.prev);
            c.next = solution.get_next_vertex(c.v);

            const auto c_v_next = this->instance.get_cost(c.v, c.next);
            c.seqrem = -this->instance.get_cost(c.prevprev, c.prev) - c_v_next;

            return c;
        }

        inline Cache1 prepare_cache1(const cobra::Solution &solution, int vertex, int backup) {
            auto c = Cache1();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v, backup);
            c.prev = solution.get_prev_vertex(route, c.v);
            c.prevprev = solution.get_prev_vertex(route, c.prev);
            c.next = solution.get_next_vertex(route, c.v);

            const auto c_v_next = this->instance.get_cost(c.v, c.next);
            c.seqrem = -this->instance.get_cost(c.prevprev, c.prev) - c_v_next;

            return c;
        }

        struct Cache2 {
            int v, next, nextnext;
            float prevrem;
        };

        inline Cache2 prepare_cache2(const cobra::Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache2();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v);
            c.next = solution.get_next_vertex(c.v);
            c.nextnext = solution.get_next_vertex(route, c.next);

            const auto c_v_next = this->instance.get_cost(c.v, c.next);
            c.prevrem = -c_v_next - this->instance.get_cost(c.next, c.nextnext);

            return c;
        }

        inline Cache2 prepare_cache2(const cobra::Solution &solution, int vertex, int backup) {
            auto c = Cache2();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v, backup);
            c.next = solution.get_next_vertex(route, c.v);
            c.nextnext = solution.get_next_vertex(route, c.next);

            const auto c_v_next = this->instance.get_cost(c.v, c.next);
            c.prevrem = -c_v_next - this->instance.get_cost(c.next, c.nextnext);
            
            return c;
        }

        template <typename C1, typename C2>
        inline float compute_cost(const C1 i, const C2 j) {

            const auto iSequenceAdd = +this->instance.get_cost(j.nextnext, i.prev) + this->instance.get_cost(i.v, j.v);
            const auto jNextAdd = +this->instance.get_cost(i.prevprev, j.next) + this->instance.get_cost(j.next, i.next);

            const auto delta = iSequenceAdd + jNextAdd + i.seqrem + j.prevrem;

            return delta;
        }
    };

}  // namespace cobra

#endif