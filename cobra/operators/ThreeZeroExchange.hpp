#ifndef _F4D_THREEZEROEXCHANGE_HPP_
#define _F4D_THREEZEROEXCHANGE_HPP_

#include "AbstractOperator.hpp"

namespace cobra {

    template <bool handle_partial_solutions = false>
    class ThreeZeroExchange : public AbstractOperator {
    public:
        ThreeZeroExchange(const Instance &instance_, MoveGenerators &moves_, float tolerance_) : AbstractOperator(instance_, moves_, tolerance_) { }

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
            const auto iPrevPrevPrev = solution.get_prev_vertex(iRoute, iPrevPrev);

            const auto jPrev = solution.get_prev_vertex(jRoute, j);

            const auto iSequenceRem = -this->instance.get_cost(iPrevPrevPrev, iPrevPrev) - this->instance.get_cost(i, iNext);
            const auto jSequenceRem = -this->instance.get_cost(j, jPrev);
            const auto iSequenceAdd = +this->instance.get_cost(jPrev, iPrevPrev) + this->instance.get_cost(i, j);
            const auto iFilling = +this->instance.get_cost(iPrevPrevPrev, iNext);

            return iSequenceAdd + iFilling + iSequenceRem + jSequenceRem;
        }

        bool is_feasible(const Solution &solution, const MoveGenerator &move) override {
            const auto i = move.get_first_vertex();
            const auto j = move.get_second_vertex();

            const auto iRoute = solution.get_route_index(i, j);
            const auto jRoute = solution.get_route_index(j, i);

            const auto iPrev = solution.get_prev_vertex(iRoute, i);
            const auto iPrevPrev = solution.get_prev_vertex(iRoute, iPrev);


            return (iRoute != jRoute && iPrev != this->instance.get_depot() && iPrevPrev != this->instance.get_depot() &&
                    solution.get_route_load(jRoute) + this->instance.get_demand(i) + this->instance.get_demand(iPrev) + this->instance.get_demand(iPrevPrev) <=
                        this->instance.get_vehicle_capacity()) ||
                   (iRoute == jRoute && j != iPrev && j != iPrevPrev && j != solution.get_next_vertex(iRoute, i));
        }

        inline void execute(Solution &solution, const MoveGenerator &move, VertexSet &storage) override {
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

            const auto jPrev = solution.get_prev_vertex(jRoute, j);

            const auto jNext = solution.get_next_vertex(jRoute, j);
            const auto jNextNext = solution.get_next_vertex(jRoute, jNext);


            storage.insert(iPrevPrevPrev);
            storage.insert(iPrevPrev);
            storage.insert(iPrev);
            storage.insert(i);
            storage.insert(iNext);
            storage.insert(iNextNext);
            storage.insert(iNextNextNext);
            storage.insert(jPrev);
            storage.insert(j);
            storage.insert(jNext);
            storage.insert(jNextNext);

            this->update_bits.at(iPrevPrevPrev, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iPrevPrev, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iPrevPrev, UPDATE_BITS_SECOND, true);
            this->update_bits.at(iPrev, UPDATE_BITS_FIRST, true);
            this->update_bits.at(i, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iNext, UPDATE_BITS_SECOND, true);
            this->update_bits.at(iNextNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iNextNextNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(jPrev, UPDATE_BITS_FIRST, true);
            this->update_bits.at(j, UPDATE_BITS_FIRST, true);
            this->update_bits.at(j, UPDATE_BITS_SECOND, true);
            this->update_bits.at(jNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(jNextNext, UPDATE_BITS_FIRST, true);


            solution.remove_vertex(iRoute, i);
            solution.remove_vertex(iRoute, iPrev);
            solution.remove_vertex(iRoute, iPrevPrev);

            solution.insert_vertex_before(jRoute, j, iPrevPrev);
            solution.insert_vertex_before(jRoute, j, iPrev);
            solution.insert_vertex_before(jRoute, j, i);

            if (solution.is_route_empty(iRoute)) {
                solution.remove_route(iRoute);
            }
        }

        void post_processing(__attribute__((unused)) Solution &solution) override { }

        std::string get_additional_statistics() override {
            return std::string();
        }

        struct Cache12 {
            int v, prev, prevprev;
            float seqrem1, seqrem2;
        };

        inline Cache12 prepare_cache12(const Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache12();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v);
            c.prev = solution.get_prev_vertex(c.v);
            c.prevprev = solution.get_prev_vertex(route, c.prev);
            const auto prevprevprev = solution.get_prev_vertex(route, c.prevprev);
            const auto next = solution.get_next_vertex(c.v);
            // filling
            c.seqrem1 = -this->instance.get_cost(prevprevprev, c.prevprev) - this->instance.get_cost(c.v, next) + this->instance.get_cost(prevprevprev, next);
            c.seqrem2 = -this->instance.get_cost(c.v, c.prev);

            return c;
        }

        inline Cache12 prepare_cache12(const Solution &solution, int vertex, int backup) {

            auto c = Cache12();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v, backup);
            c.prev = solution.get_prev_vertex(route, c.v);
            c.prevprev = solution.get_prev_vertex(route, c.prev);
            const auto prevprevprev = solution.get_prev_vertex(route, c.prevprev);
            const auto next = solution.get_next_vertex(route, c.v);
            // filling
            c.seqrem1 = -this->instance.get_cost(prevprevprev, c.prevprev) - this->instance.get_cost(c.v, next) + this->instance.get_cost(prevprevprev, next);
            c.seqrem2 = -this->instance.get_cost(c.v, c.prev);


            return c;
        }

        inline std::pair<float, float> compute_cost_pair(const struct Cache12 i, const struct Cache12 j) {

            const auto c_iv_jv = this->instance.get_cost(i.v, j.v);

            const auto iSequenceAdd = +this->instance.get_cost(j.prev, i.prevprev) + c_iv_jv;
            const auto jSequenceAdd = +this->instance.get_cost(i.prev, j.prevprev) + c_iv_jv;

            const auto delta1 = iSequenceAdd + i.seqrem1 + j.seqrem2;
            const auto delta2 = jSequenceAdd + j.seqrem1 + i.seqrem2;

            return {delta1, delta2};
        }

        struct Cache1 {
            int v, prevprev;
            float seqrem1;
        };

        inline Cache1 prepare_cache1(const Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache1();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v);
            const auto prev = solution.get_prev_vertex(c.v);
            c.prevprev = solution.get_prev_vertex(route, prev);
            const auto prevprevprev = solution.get_prev_vertex(route, c.prevprev);
            const auto next = solution.get_next_vertex(c.v);

            c.seqrem1 = -this->instance.get_cost(prevprevprev, c.prevprev) - this->instance.get_cost(c.v, next) + this->instance.get_cost(prevprevprev, next);


            return c;
        }

        inline Cache1 prepare_cache1(const Solution &solution, int vertex, int backup) {
            auto c = Cache1();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v, backup);
            const auto prev = solution.get_prev_vertex(route, c.v);
            c.prevprev = solution.get_prev_vertex(route, prev);
            const auto prevprevprev = solution.get_prev_vertex(route, c.prevprev);
            const auto next = solution.get_next_vertex(route, c.v);

            c.seqrem1 = -this->instance.get_cost(prevprevprev, c.prevprev) - this->instance.get_cost(c.v, next) + this->instance.get_cost(prevprevprev, next);


            return c;
        }

        struct Cache2 {
            int v, prev;
            float seqrem2;
        };

        inline Cache2 prepare_cache2(const Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache2();
            c.v = vertex;
            c.prev = solution.get_prev_vertex(c.v);

            c.seqrem2 = -this->instance.get_cost(c.v, c.prev);


            return c;
        }

        inline Cache2 prepare_cache2(const Solution &solution, int vertex, int backup) {
            auto c = Cache2();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v, backup);
            c.prev = solution.get_prev_vertex(route, c.v);

            c.seqrem2 = -this->instance.get_cost(c.v, c.prev);


            return c;
        }

        template <typename C1, typename C2>
        inline float compute_cost(const C1 i, const C2 j) {

            const auto iSequenceAdd = +this->instance.get_cost(j.prev, i.prevprev) + this->instance.get_cost(i.v, j.v);

            const auto delta = iSequenceAdd + i.seqrem1 + j.seqrem2;

            return delta;
        }
    };


}  // namespace cobra

#endif