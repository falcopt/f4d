#ifndef _F4D_TWOTWOEXCHANGE_HPP_
#define _F4D_TWOTWOEXCHANGE_HPP_

#include "AbstractOperator.hpp"

namespace cobra {

    template <bool handle_partial_solutions = false>
    class TwoTwoExchange : public AbstractOperator {
    public:
        TwoTwoExchange(const Instance &instance_, MoveGenerators &moves_, float tolerance_) : AbstractOperator(instance_, moves_, tolerance_) { }

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

            const auto jPrev = solution.get_prev_vertex(jRoute, j);
            const auto jPrevPrev = solution.get_prev_vertex(jRoute, jPrev);
            const auto jPrevPrevPrev = solution.get_prev_vertex(jRoute, jPrevPrev);

            const auto iSequenceRem = -this->instance.get_cost(iPrevPrev, iPrev) - this->instance.get_cost(i, iNext);
            const auto jSequenceRem = -this->instance.get_cost(jPrevPrevPrev, jPrevPrev) - this->instance.get_cost(jPrev, j);

            const auto iSequenceAdd = +this->instance.get_cost(jPrevPrevPrev, iPrev) + this->instance.get_cost(i, j);
            const auto jSequenceAdd = +this->instance.get_cost(iPrevPrev, jPrevPrev) + this->instance.get_cost(jPrev, iNext);

            const auto delta = iSequenceAdd + jSequenceAdd + iSequenceRem + jSequenceRem;

            return delta;
        }

        bool is_feasible(const Solution &solution, const MoveGenerator &move) override {
            const auto i = move.get_first_vertex();
            const auto j = move.get_second_vertex();

            const auto iRoute = solution.get_route_index(i, j);
            const auto jRoute = solution.get_route_index(j, i);

            const auto iPrev = solution.get_prev_vertex(iRoute, i);
            const auto jPrev = solution.get_prev_vertex(jRoute, j);
            const auto jPrevPrev = solution.get_prev_vertex(jRoute, jPrev);


            return (iRoute != jRoute && iPrev != this->instance.get_depot() && jPrev != this->instance.get_depot() && jPrevPrev != this->instance.get_depot() &&
                    solution.get_route_load(jRoute) - this->instance.get_demand(jPrev) - this->instance.get_demand(jPrevPrev) + this->instance.get_demand(i) +
                            this->instance.get_demand(iPrev) <=
                        this->instance.get_vehicle_capacity() &&
                    solution.get_route_load(iRoute) + this->instance.get_demand(jPrev) + this->instance.get_demand(jPrevPrev) - this->instance.get_demand(i) -
                            this->instance.get_demand(iPrev) <=
                        this->instance.get_vehicle_capacity()) ||
                   (iRoute == jRoute && i != jPrev && i != jPrevPrev && solution.get_next_vertex(iRoute, i) != jPrevPrev && j != iPrev);
        }


        inline void execute(Solution &solution, const MoveGenerator &move, VertexSet &storage) override {

            const auto i = move.get_first_vertex();
            const auto j = move.get_second_vertex();

            const auto iRoute = solution.get_route_index(i, j);
            const auto jRoute = solution.get_route_index(j, i);

            const auto iPrev = solution.get_prev_vertex(iRoute, i);
            const auto iNext = solution.get_next_vertex(iRoute, i);
            const auto iPrevPrev = solution.get_prev_vertex(iRoute, iPrev);

            const auto iNextNext = solution.get_next_vertex(iRoute, iNext);
            const auto iNextNextNext = solution.get_next_vertex(iRoute, iNextNext);

            const auto jPrev = solution.get_prev_vertex(jRoute, j);
            const auto jPrevPrev = solution.get_prev_vertex(jRoute, jPrev);
            const auto jPrevPrevPrev = solution.get_prev_vertex(jRoute, jPrevPrev);

            const auto jNext = solution.get_next_vertex(jRoute, j);
            const auto jNextNext = solution.get_next_vertex(jRoute, jNext);


            storage.insert(iPrevPrev);
            storage.insert(iPrev);
            storage.insert(i);
            storage.insert(iNext);
            storage.insert(iNextNext);
            storage.insert(iNextNextNext);
            storage.insert(jPrevPrevPrev);
            storage.insert(jPrevPrev);
            storage.insert(jPrev);
            storage.insert(j);
            storage.insert(jNext);
            storage.insert(jNextNext);

            this->update_bits.at(iPrevPrev, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iPrev, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iPrev, UPDATE_BITS_SECOND, true);
            this->update_bits.at(i, UPDATE_BITS_FIRST, true);
            this->update_bits.at(i, UPDATE_BITS_SECOND, true);
            this->update_bits.at(iNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iNext, UPDATE_BITS_SECOND, true);
            this->update_bits.at(iNextNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iNextNext, UPDATE_BITS_SECOND, true);
            this->update_bits.at(iNextNextNext, UPDATE_BITS_SECOND, true);
            this->update_bits.at(jPrevPrevPrev, UPDATE_BITS_FIRST, true);
            this->update_bits.at(jPrevPrev, UPDATE_BITS_FIRST, true);
            this->update_bits.at(jPrevPrev, UPDATE_BITS_SECOND, true);
            this->update_bits.at(jPrev, UPDATE_BITS_FIRST, true);
            this->update_bits.at(jPrev, UPDATE_BITS_SECOND, true);
            this->update_bits.at(j, UPDATE_BITS_FIRST, true);
            this->update_bits.at(j, UPDATE_BITS_SECOND, true);
            this->update_bits.at(jNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(jNext, UPDATE_BITS_SECOND, true);
            this->update_bits.at(jNextNext, UPDATE_BITS_SECOND, true);


            solution.remove_vertex(iRoute, i);
            solution.remove_vertex(iRoute, iPrev);

            solution.insert_vertex_before(jRoute, j, iPrev);
            solution.insert_vertex_before(jRoute, j, i);

            solution.remove_vertex(jRoute, jPrev);
            solution.remove_vertex(jRoute, jPrevPrev);
            solution.insert_vertex_before(iRoute, iNext, jPrevPrev);
            solution.insert_vertex_before(iRoute, iNext, jPrev);
        }

        void post_processing(__attribute__((unused)) Solution &solution) override { }

        std::string get_additional_statistics() override {
            return std::string();
        }

        struct Cache12 {
            int v, prev, prevprev, prevprevprev, next;
            float seqrem1, seqrem2;
        };

        inline Cache12 prepare_cache12(const Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache12();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v);
            c.prev = solution.get_prev_vertex(c.v);
            c.prevprev = solution.get_prev_vertex(route, c.prev);
            c.prevprevprev = solution.get_prev_vertex(route, c.prevprev);
            c.next = solution.get_next_vertex(c.v);

            c.seqrem1 = -this->instance.get_cost(c.prevprev, c.prev) - this->instance.get_cost(c.v, c.next);
            c.seqrem2 = -this->instance.get_cost(c.prevprevprev, c.prevprev) - this->instance.get_cost(c.prev, c.v);


            return c;
        }


        inline Cache12 prepare_cache12(const Solution &solution, int vertex, int backup) {

            auto c = Cache12();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v, backup);
            c.prev = solution.get_prev_vertex(route, c.v);
            c.prevprev = solution.get_prev_vertex(route, c.prev);
            c.prevprevprev = solution.get_prev_vertex(route, c.prevprev);
            c.next = solution.get_next_vertex(route, c.v);

            c.seqrem1 = -this->instance.get_cost(c.prevprev, c.prev) - this->instance.get_cost(c.v, c.next);
            c.seqrem2 = -this->instance.get_cost(c.prevprevprev, c.prevprev) - this->instance.get_cost(c.prev, c.v);


            return c;
        }

        inline std::pair<float, float> compute_cost_pair(const struct Cache12 i, const struct Cache12 j) {

            const auto c_iv_jv = this->instance.get_cost(i.v, j.v);
            const auto c_iprevprev_jprevprev = this->instance.get_cost(i.prevprev, j.prevprev);

            const auto iSequenceAdd1 = +this->instance.get_cost(j.prevprevprev, i.prev) + c_iv_jv;
            const auto jSequenceAdd2 = +this->instance.get_cost(i.prevprevprev, j.prev) + c_iv_jv;

            const auto jSequenceAdd1 = +c_iprevprev_jprevprev + this->instance.get_cost(j.prev, i.next);
            const auto iSequenceAdd2 = +c_iprevprev_jprevprev + this->instance.get_cost(i.prev, j.next);

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

            c.seqrem1 = -this->instance.get_cost(c.prevprev, c.prev) - this->instance.get_cost(c.v, c.next);


            return c;
        }

        inline Cache1 prepare_cache1(const Solution &solution, int vertex, int backup) {
            auto c = Cache1();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v, backup);
            c.prev = solution.get_prev_vertex(route, c.v);
            c.prevprev = solution.get_prev_vertex(route, c.prev);
            c.next = solution.get_next_vertex(route, c.v);

            c.seqrem1 = -this->instance.get_cost(c.prevprev, c.prev) - this->instance.get_cost(c.v, c.next);


            return c;
        }

        struct Cache2 {
            int v, prev, prevprev, prevprevprev;
            float seqrem2;
        };

        inline Cache2 prepare_cache2(const Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache2();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v);
            c.prev = solution.get_prev_vertex(c.v);
            c.prevprev = solution.get_prev_vertex(route, c.prev);
            c.prevprevprev = solution.get_prev_vertex(route, c.prevprev);

            c.seqrem2 = -this->instance.get_cost(c.prevprevprev, c.prevprev) - this->instance.get_cost(c.prev, c.v);


            return c;
        }

        inline Cache2 prepare_cache2(const Solution &solution, int vertex, int backup) {
            auto c = Cache2();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v, backup);
            c.prev = solution.get_prev_vertex(route, c.v);
            c.prevprev = solution.get_prev_vertex(route, c.prev);
            c.prevprevprev = solution.get_prev_vertex(route, c.prevprev);

            c.seqrem2 = -this->instance.get_cost(c.prevprevprev, c.prevprev) - this->instance.get_cost(c.prev, c.v);


            return c;
        }

        template <typename C1, typename C2>
        inline float compute_cost(const C1 i, const C2 j) {

            const auto iSequenceAdd = +this->instance.get_cost(j.prevprevprev, i.prev) + this->instance.get_cost(i.v, j.v);
            const auto jSequenceAdd = +this->instance.get_cost(i.prevprev, j.prevprev) + this->instance.get_cost(j.prev, i.next);

            const auto delta = iSequenceAdd + jSequenceAdd + i.seqrem1 + j.seqrem2;

            return delta;
        }
    };

}  // namespace cobra

#endif  