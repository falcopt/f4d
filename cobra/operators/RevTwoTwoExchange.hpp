#ifndef _F4D_REVTWOTWOEXCHANGE_HPP_
#define _F4D_REVTWOTWOEXCHANGE_HPP_

#include "AbstractOperator.hpp"

namespace cobra {

    template <bool handle_partial_solutions = false, bool reverse_both_strings = false>
    class RevTwoTwoExchange : public AbstractOperator {
    public:
        RevTwoTwoExchange(const Instance &instance_, MoveGenerators &moves_, float tolerance_) : AbstractOperator(instance_, moves_, tolerance_) { }

        static constexpr bool is_symmetric = false;

    protected:
        inline void pre_processing(__attribute__((unused)) Solution &solution) override { }

        inline float compute_cost(const Solution &solution, const MoveGenerator &move) override {
            const auto i = move.get_first_vertex();
            const auto j = move.get_second_vertex();

            const auto iRoute = solution.get_route_index(i, j);
            const auto jRoute = solution.get_route_index(j, i);

            const auto iNext = solution.get_next_vertex(iRoute, i);
            const auto iPrev = solution.get_prev_vertex(iRoute, i);
            const auto iPrevPrev = solution.get_prev_vertex(iRoute, iPrev);

            const auto jNext = solution.get_next_vertex(jRoute, j);
            const auto jNextNext = solution.get_next_vertex(jRoute, jNext);
            const auto jNextNextNext = solution.get_next_vertex(jRoute, jNextNext);

            const auto iSequenceRem = -solution.get_prev_cost(iPrev) - solution.get_next_cost(i);
            const auto jSequenceRem = -solution.get_next_cost(j) - solution.get_next_cost(jNextNext);

            const auto iSequenceAdd = +this->instance.get_cost(jNextNextNext, iPrev) + this->instance.get_cost(i, j);

            float jSequenceAdd;
            if constexpr (reverse_both_strings) {
                jSequenceAdd = +this->instance.get_cost(iPrevPrev, jNextNext) + this->instance.get_cost(jNext, iNext);
            } else {
                jSequenceAdd = +this->instance.get_cost(iPrevPrev, jNext) + this->instance.get_cost(jNextNext, iNext);
            }

            const auto delta = iSequenceAdd + jSequenceAdd + iSequenceRem + jSequenceRem;

            return delta;
        }

        bool is_feasible(const Solution &solution, const MoveGenerator &move) override {
            const auto i = move.get_first_vertex();
            const auto j = move.get_second_vertex();

            const auto iRoute = solution.get_route_index(i, j);
            const auto jRoute = solution.get_route_index(j, i);

            const auto iPrev = solution.get_prev_vertex(iRoute, i);
            const auto jNext = solution.get_next_vertex(jRoute, j);
            const auto jNextNext = solution.get_next_vertex(jRoute, jNext);


            return (iRoute != jRoute && iPrev != this->instance.get_depot() && jNext != this->instance.get_depot() && jNextNext != this->instance.get_depot() &&
                    solution.get_route_load(jRoute) - this->instance.get_demand(jNext) - this->instance.get_demand(jNextNext) + this->instance.get_demand(i) +
                            this->instance.get_demand(iPrev) <=
                        this->instance.get_vehicle_capacity() &&
                    solution.get_route_load(iRoute) + this->instance.get_demand(jNext) + this->instance.get_demand(jNextNext) - this->instance.get_demand(i) -
                            this->instance.get_demand(iPrev) <=
                        this->instance.get_vehicle_capacity()) ||
                   (iRoute == jRoute && j != iPrev && jNext != iPrev && jNextNext != iPrev && jNextNext != solution.get_prev_vertex(iRoute, iPrev));
        }

        inline void execute(Solution &solution, const MoveGenerator &move, VertexSet &storage) override {

            const auto i = move.get_first_vertex();
            const auto j = move.get_second_vertex();

            const auto iRoute = solution.get_route_index(i, j);
            const auto jRoute = solution.get_route_index(j, i);

            const auto iPrev = solution.get_prev_vertex(iRoute, i);
            const auto iNext = solution.get_next_vertex(iRoute, i);
            const auto iPrevPrev = solution.get_prev_vertex(iRoute, iPrev);
            const auto iPrevPrevPrev = solution.get_prev_vertex(iRoute, iPrevPrev);
            const auto iPrevPrevPrevPrev = solution.get_prev_vertex(iRoute, iPrevPrevPrev);

            const auto iNextNext = solution.get_next_vertex(iRoute, iNext);

            const auto jPrev = solution.get_prev_vertex(jRoute, j);
            const auto jPrevPrev = solution.get_prev_vertex(jRoute, jPrev);
            const auto jNext = solution.get_next_vertex(jRoute, j);
            const auto jNextNext = solution.get_next_vertex(jRoute, jNext);
            const auto jNextNextNext = solution.get_next_vertex(jRoute, jNextNext);
            const auto jNextNextNextNext = solution.get_next_vertex(jRoute, jNextNextNext);


            storage.insert(iPrevPrevPrevPrev);
            storage.insert(iPrevPrevPrev);
            storage.insert(iPrevPrev);
            storage.insert(iPrev);
            storage.insert(i);
            storage.insert(iNext);
            storage.insert(iNextNext);
            storage.insert(jPrevPrev);
            storage.insert(jPrev);
            storage.insert(j);
            storage.insert(jNext);
            storage.insert(jNextNext);
            storage.insert(jNextNextNext);
            storage.insert(jNextNextNextNext);

            this->update_bits.at(iPrevPrevPrevPrev, UPDATE_BITS_SECOND, true);
            this->update_bits.at(iPrevPrevPrev, UPDATE_BITS_SECOND, true);
            this->update_bits.at(iPrevPrev, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iPrevPrev, UPDATE_BITS_SECOND, true);
            this->update_bits.at(iPrev, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iPrev, UPDATE_BITS_SECOND, true);
            this->update_bits.at(i, UPDATE_BITS_FIRST, true);
            this->update_bits.at(i, UPDATE_BITS_SECOND, true);
            this->update_bits.at(iNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iNextNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(jNextNextNextNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(jNextNextNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(jNextNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(jNextNext, UPDATE_BITS_SECOND, true);
            this->update_bits.at(jNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(jNext, UPDATE_BITS_SECOND, true);
            this->update_bits.at(j, UPDATE_BITS_FIRST, true);
            this->update_bits.at(j, UPDATE_BITS_SECOND, true);
            this->update_bits.at(jPrev, UPDATE_BITS_SECOND, true);
            this->update_bits.at(jPrevPrev, UPDATE_BITS_SECOND, true);


            solution.remove_vertex(iRoute, i);
            solution.remove_vertex(iRoute, iPrev);

            solution.insert_vertex_before(jRoute, jNextNextNext, i);
            solution.insert_vertex_before(jRoute, jNextNextNext, iPrev);

            solution.remove_vertex(jRoute, jNext);
            solution.remove_vertex(jRoute, jNextNext);

            if constexpr (reverse_both_strings) {
                solution.insert_vertex_before(iRoute, iNext, jNextNext);
                solution.insert_vertex_before(iRoute, iNext, jNext);
            } else {
                solution.insert_vertex_before(iRoute, iNext, jNext);
                solution.insert_vertex_before(iRoute, iNext, jNextNext);
            }
        }

        void post_processing(__attribute__((unused)) Solution &solution) override { }

        std::string get_additional_statistics() override {
            return std::string();
        }

        struct Cache12 {
            int v, prev, prevprev, next, nextnext, nextnextnext;
            float seqrem1, seqrem2;
        };

        inline Cache12 prepare_cache12(const Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache12();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v);
            c.prev = solution.get_prev_vertex(c.v);
            c.prevprev = solution.get_prev_vertex(route, c.prev);
            c.next = solution.get_next_vertex(c.v);
            c.nextnext = solution.get_next_vertex(route, c.next);
            c.nextnextnext = solution.get_next_vertex(route, c.nextnext);

            const auto c_v_next = solution.get_next_cost(c.v);
            c.seqrem1 = -solution.get_prev_cost(c.prev, c.prevprev) - c_v_next;
            c.seqrem2 = -c_v_next - solution.get_next_cost(c.nextnext, c.nextnextnext);


            return c;
        }

        inline Cache12 prepare_cache12(const Solution &solution, int vertex, int backup) {

            auto c = Cache12();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v, backup);
            c.prev = solution.get_prev_vertex(route, c.v);
            c.prevprev = solution.get_prev_vertex(route, c.prev);
            c.next = solution.get_next_vertex(route, c.v);
            c.nextnext = solution.get_next_vertex(route, c.next);
            c.nextnextnext = solution.get_next_vertex(route, c.nextnext);

            const auto c_v_next = solution.get_next_cost(c.v, c.next);
            c.seqrem1 = -solution.get_prev_cost(c.prev, c.prevprev) - c_v_next;
            c.seqrem2 = -c_v_next - solution.get_next_cost(c.nextnext, c.nextnextnext);


            return c;
        }

        inline std::pair<float, float> compute_cost_pair(const struct Cache12 i, const struct Cache12 j) {

            const auto c_iv_jv = this->instance.get_cost(i.v, j.v);

            const auto iSequenceAdd1 = +this->instance.get_cost(j.nextnextnext, i.prev) + c_iv_jv;
            const auto jSequenceAdd2 = +this->instance.get_cost(i.nextnextnext, j.prev) + c_iv_jv;

            float jSequenceAdd1, iSequenceAdd2;
            if constexpr (reverse_both_strings) {
                const auto c_inext_jnext = this->instance.get_cost(i.next, j.next);
                jSequenceAdd1 = +this->instance.get_cost(i.prevprev, j.nextnext) + c_inext_jnext;
                iSequenceAdd2 = +this->instance.get_cost(j.prevprev, i.nextnext) + c_inext_jnext;
            } else {
                jSequenceAdd1 = +this->instance.get_cost(i.prevprev, j.next) + this->instance.get_cost(j.nextnext, i.next);
                iSequenceAdd2 = +this->instance.get_cost(j.prevprev, i.next) + this->instance.get_cost(i.nextnext, j.next);
            }


            const auto delta1 = iSequenceAdd1 + jSequenceAdd1 + i.seqrem1 + j.seqrem2;
            const auto delta2 = jSequenceAdd2 + iSequenceAdd2 + j.seqrem1 + i.seqrem2;

            return {delta1, delta2};
        }

        struct Cache1 {
            int v, prev, prevprev, next;
            float seqrem1;
        };

        inline Cache1 prepare_cache1(const Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache1();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v);
            c.prev = solution.get_prev_vertex(c.v);
            c.prevprev = solution.get_prev_vertex(route, c.prev);
            c.next = solution.get_next_vertex(c.v);

            const auto c_v_next = solution.get_next_cost(c.v);
            c.seqrem1 = -solution.get_prev_cost(c.prev, c.prevprev) - c_v_next;


            return c;
        }

        inline Cache1 prepare_cache1(const Solution &solution, int vertex, int backup) {
            auto c = Cache1();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v, backup);
            c.prev = solution.get_prev_vertex(route, c.v);
            c.prevprev = solution.get_prev_vertex(route, c.prev);
            c.next = solution.get_next_vertex(route, c.v);

            const auto c_v_next = solution.get_next_cost(c.v, c.next);
            c.seqrem1 = -solution.get_prev_cost(c.prev, c.prevprev) - c_v_next;


            return c;
        }

        struct Cache2 {
            int v, next, nextnext, nextnextnext;
            float seqrem2;
        };

        inline Cache2 prepare_cache2(const Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache2();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v);
            c.next = solution.get_next_vertex(c.v);
            c.nextnext = solution.get_next_vertex(route, c.next);
            c.nextnextnext = solution.get_next_vertex(route, c.nextnext);

            const auto c_v_next = solution.get_next_cost(c.v);
            c.seqrem2 = -c_v_next - solution.get_next_cost(c.nextnext, c.nextnextnext);


            return c;
        }

        inline Cache2 prepare_cache2(const Solution &solution, int vertex, int backup) {
            auto c = Cache2();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v, backup);
            c.next = solution.get_next_vertex(route, c.v);
            c.nextnext = solution.get_next_vertex(route, c.next);
            c.nextnextnext = solution.get_next_vertex(route, c.nextnext);

            const auto c_v_next = solution.get_next_cost(c.v, c.next);
            c.seqrem2 = -c_v_next - solution.get_next_cost(c.nextnext, c.nextnextnext);


            return c;
        }

        template <typename C1, typename C2>
        inline float compute_cost(const C1 i, const C2 j) {

            const auto iSequenceAdd = +this->instance.get_cost(j.nextnextnext, i.prev) + this->instance.get_cost(i.v, j.v);

            float jSequenceAdd;
            if constexpr (reverse_both_strings) {
                jSequenceAdd = +this->instance.get_cost(i.prevprev, j.nextnext) + this->instance.get_cost(j.next, i.next);
            } else {
                jSequenceAdd = +this->instance.get_cost(i.prevprev, j.next) + this->instance.get_cost(j.nextnext, i.next);
            }

            const auto delta = iSequenceAdd + jSequenceAdd + i.seqrem1 + j.seqrem2;

            return delta;
        }
    };

}  // namespace cobra

#endif