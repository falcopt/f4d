#ifndef _F4D_SPLITEXCHANGE_HPP_
#define _F4D_SPLITEXCHANGE_HPP_

#include "AbstractOperator.hpp"

namespace cobra {

    template <bool handle_partial_solutions = false>
    class SplitExchange : public AbstractOperator {
    public:
        SplitExchange(const Instance &instance_, MoveGenerators &moves_, float tolerance_) : AbstractOperator(instance_, moves_, tolerance_) { }

        static constexpr bool is_symmetric = true;

    protected:
        inline void pre_processing(Solution &solution) override {

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
            const auto jNext = solution.get_next_vertex(jRoute, j);

            return -this->instance.get_cost(i, iNext) + this->instance.get_cost(i, j) - this->instance.get_cost(j, jNext) +
                   this->instance.get_cost(jNext, iNext);
        }

        bool is_feasible(const Solution &solution, const MoveGenerator &move) override {
            const auto i = move.get_first_vertex();
            const auto j = move.get_second_vertex();

            const auto iRoute = solution.get_route_index(i, j);
            const auto jRoute = solution.get_route_index(j, i);


            return iRoute != jRoute &&
                   (solution.get_route_load_before_included(i) + solution.get_route_load_before_included(j) <= this->instance.get_vehicle_capacity() &&
                    solution.get_route_load_after_included(j) - this->instance.get_demand(j) + solution.get_route_load_after_included(i) -
                            this->instance.get_demand(i) <=
                        this->instance.get_vehicle_capacity());
        }

        inline void execute(Solution &solution, const MoveGenerator &move, VertexSet &storage) override {

            const auto i = move.get_first_vertex();
            const auto j = move.get_second_vertex();

            const auto iRoute = solution.get_route_index(i, j);
            const auto jRoute = solution.get_route_index(j, i);

            assert(solution.get_first_customer(iRoute) != this->instance.get_depot());
            assert(solution.get_first_customer(jRoute) != this->instance.get_depot());

            storage.insert(this->instance.get_depot());
            for (auto curr = i; curr != this->instance.get_depot(); curr = solution.get_next_vertex(curr)) {
                storage.insert(curr);
            }

            const auto jNextNext = solution.get_next_vertex(jRoute, solution.get_next_vertex(j));
            const auto jStop = jNextNext == solution.get_first_customer(jRoute) ? this->instance.get_depot() : jNextNext;  // handle special case
            for (auto curr = solution.get_first_customer(jRoute); curr != jStop; curr = solution.get_next_vertex(curr)) {
                storage.insert(curr);
            }


            solution.split(i, iRoute, j, jRoute);

            if (solution.is_route_empty(iRoute)) {
                solution.remove_route(iRoute);
            } else {
                solution.update_cumulative_route_loads(iRoute);
            }

            if (solution.is_route_empty(jRoute)) {
                solution.remove_route(jRoute);
            } else {
                solution.update_cumulative_route_loads(jRoute);
            }
        }

        void post_processing(__attribute__((unused)) Solution &solution) override { }

        std::string get_additional_statistics() override {
            return std::string();
        }

        struct Cache12 {
            int v, next;
            float seqrem;
        };

        inline Cache12 prepare_cache12(const Solution &solution, int vertex) {
            assert(vertex != this->instance.get_depot());
            auto c = Cache12();
            c.v = vertex;
            c.next = solution.get_next_vertex(c.v);

            c.seqrem = -this->instance.get_cost(c.v, c.next);

            return c;
        }

        inline Cache12 prepare_cache12(const Solution &solution, int vertex, int backup) {

            auto c = Cache12();
            c.v = vertex;
            const auto route = solution.get_route_index(c.v, backup);
            c.next = solution.get_next_vertex(route, c.v);

            c.seqrem = -this->instance.get_cost(c.v, c.next);

            return c;
        }


        inline float compute_cost(const struct Cache12 i, const struct Cache12 j) {

            const auto iSequenceAdd = this->instance.get_cost(i.v, j.v) + this->instance.get_cost(j.next, i.next);

            const auto delta = iSequenceAdd + i.seqrem + j.seqrem;

            return delta;
        }
    };

}  // namespace cobra

#endif