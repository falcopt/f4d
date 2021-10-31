#ifndef _F4D_ONEZEROEXCHANGE_HPP_
#define _F4D_ONEZEROEXCHANGE_HPP_

#include "AbstractOperator.hpp"

namespace cobra {

    template <bool handle_partial_solutions = false>
    class OneZeroExchange : public AbstractOperator {

    public:
        OneZeroExchange(const Instance &instance_, MoveGenerators &moves_, float tolerance_) : AbstractOperator(instance_, moves_, tolerance_) { }

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

            return -this->instance.get_cost(iPrev, i) - this->instance.get_cost(i, iNext) + this->instance.get_cost(iPrev, iNext) -
                   this->instance.get_cost(jPrev, j) + this->instance.get_cost(jPrev, i) + this->instance.get_cost(i, j);
        }

        bool is_feasible(const Solution &solution, const MoveGenerator &move) override {
            const auto i = move.get_first_vertex();
            const auto j = move.get_second_vertex();

            const auto iRoute = solution.get_route_index(i, j);
            const auto jRoute = solution.get_route_index(j, i);

            return (iRoute != jRoute && solution.get_route_load(jRoute) + this->instance.get_demand(i) <= this->instance.get_vehicle_capacity()) ||
                   (iRoute == jRoute && j != solution.get_next_vertex(iRoute, i));
        }

        inline void execute(Solution &solution, const MoveGenerator &move, VertexSet &storage) override {

            const auto i = move.get_first_vertex();
            const auto j = move.get_second_vertex();

            const auto iRoute = solution.get_route_index(i, j);
            const auto jRoute = solution.get_route_index(j, i);

            const auto iPrev = solution.get_prev_vertex(iRoute, i);
            const auto iNext = solution.get_next_vertex(iRoute, i);

            const auto jPrev = solution.get_prev_vertex(jRoute, j);

            storage.insert(iPrev);
            storage.insert(i);
            storage.insert(iNext);
            storage.insert(jPrev);
            storage.insert(j);

            this->update_bits.at(iPrev, UPDATE_BITS_FIRST, true);
            this->update_bits.at(i, UPDATE_BITS_FIRST, true);
            this->update_bits.at(i, UPDATE_BITS_SECOND, true);
            this->update_bits.at(iNext, UPDATE_BITS_FIRST, true);
            this->update_bits.at(iNext, UPDATE_BITS_SECOND, true);
            this->update_bits.at(j, UPDATE_BITS_FIRST, true);
            this->update_bits.at(j, UPDATE_BITS_SECOND, true);
            this->update_bits.at(jPrev, UPDATE_BITS_FIRST, true);


            solution.remove_vertex(iRoute, i);
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

            const auto delta1 = i.vrem + j.prevrem + this->instance.get_cost(j.prev, i.v) + c_iv_jv;
            const auto delta2 = j.vrem + i.prevrem + this->instance.get_cost(i.prev, j.v) + c_iv_jv;

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
            return i.vrem + j.prevrem + this->instance.get_cost(j.prev, i.v) + this->instance.get_cost(i.v, j.v);
        }
    };

}  // namespace cobra

#endif