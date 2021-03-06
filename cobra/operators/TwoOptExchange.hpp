#ifndef _F4D_TWOOPTEXCHANGE_HPP_
#define _F4D_TWOOPTEXCHANGE_HPP_

#include "AbstractOperator.hpp"

namespace cobra {

    template <bool handle_partial_solutions = false>
    class TwoOptExchange : public AbstractOperator {
    public:
        TwoOptExchange(const Instance &instance_, MoveGenerators &moves_, float tolerance_) : AbstractOperator(instance_, moves_, tolerance_) { }

        static constexpr bool is_symmetric = true;

    protected:
        inline void pre_processing(__attribute__((unused)) Solution &solution) override { }

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

            return iRoute == jRoute;
        }

        inline void execute(Solution &solution, const MoveGenerator &move, VertexSet &storage) override {

            const auto i = move.get_first_vertex();
            const auto j = move.get_second_vertex();

            const auto iRoute = solution.get_route_index(i, j);

            assert(solution.get_first_customer(iRoute) != this->instance.get_depot());

            const auto jNextNext = solution.get_next_vertex(iRoute, solution.get_next_vertex(iRoute, j));
            // selective update of i, iNext, j, jNext that are directly involved and of the reversed part only in which prev & next pointer where changed
            // use the do-while for very short tour (4 vertices) in which jNextNext equal i
            auto curr = i;
            do {
                storage.insert(curr);
                curr = solution.get_next_vertex(iRoute, curr);
            } while (curr != jNextNext);

            const auto iNext = solution.get_next_vertex(iRoute, i);

            solution.reverse_route_path(iRoute, iNext, j);
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
            return iSequenceAdd + i.seqrem + j.seqrem;
        }
    };

}  // namespace cobra

#endif