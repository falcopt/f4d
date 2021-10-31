#ifndef _F4D_LOCALSEARCH_HPP_
#define _F4D_LOCALSEARCH_HPP_

#include <cfloat>
#include <map>
#include <set>

#include "Instance.hpp"
#include "MoveGenerators.hpp"
#include "Solution.hpp"
#include "VertexSet.hpp"
#include "operators/AbstractOperator.hpp"
#include "operators/EjectionChain.hpp"
#include "operators/OneOneExchange.hpp"
#include "operators/OneZeroExchange.hpp"
#include "operators/RevThreeOneExchange.hpp"
#include "operators/RevThreeThreeExchange.hpp"
#include "operators/RevThreeTwoExchange.hpp"
#include "operators/RevThreeZeroExchange.hpp"
#include "operators/RevTwoOneExchange.hpp"
#include "operators/RevTwoTwoExchange.hpp"
#include "operators/RevTwoZeroExchange.hpp"
#include "operators/SplitExchange.hpp"
#include "operators/TailsExchange.hpp"
#include "operators/ThreeOneExchange.hpp"
#include "operators/ThreeThreeExchange.hpp"
#include "operators/ThreeTwoExchange.hpp"
#include "operators/ThreeZeroExchange.hpp"
#include "operators/TwoOneExchange.hpp"
#include "operators/TwoOptExchange.hpp"
#include "operators/TwoTwoExchange.hpp"
#include "operators/TwoZeroExchange.hpp"
#include "operators/TailChain.hpp"
#include "operators/SplitChain.hpp"

namespace cobra {

    /**
     * Supported local search operators
     */
    enum Operator { E10, E11, E20, E21, E22, E30, E31, E32, E33, SPLIT, TAILS, TWOPT, EJCH, RE20, RE21, RE22B, RE22S, RE30, RE31, RE32B, RE32S, RE33B, RE33S, TLCH, STCH };

    class VariableNeighborhoodDescentInterface {

    public:
        virtual void apply(Solution& solution) = 0;
    };


    /**
     * Local search statistics
     */
    template <bool handle_partial_solutions, bool log_statistics>
    class RandomizedVariableNeighborhoodDescent;  // forward declaration

    class Statistics {

        template <bool, bool>
        friend class RandomizedVariableNeighborhoodDescent;

        std::vector<std::pair<float,                                     // base objective value (e.g. current solution value)
                              std::vector<std::pair<Operator, float>>>>  // list of operator id, objective value after operator application
            data;

        // temporary variables
        float base_obj{};
        std::vector<std::pair<Operator, float>> sequence;

        void begin(const cobra::Solution& solution) { base_obj = solution.get_cost(); }

        void append(Operator op, const cobra::Solution& solution) { sequence.emplace_back(op, solution.get_cost()); }

        void end() {
            data.emplace_back(base_obj, sequence);
            sequence.clear();
        }

    public:
        void list() {}
    };

    /**
     * Templates to enable/disable statistics class members.
     */
    template <bool, bool handle_partial_solutions>
    struct EnableStatistics {};
    template <bool handle_partial_solutions>
    struct EnableStatistics<true, handle_partial_solutions> {

        std::unordered_map<AbstractOperator*, Operator> operator_ptr_to_id;

        Statistics stats;
    };

    /**
     * Local search manager
     * @tparam handle_partial_solutions, when this is true, local search can be safely applied yo
     * solutions where some customers are not served.
     * @tparam log_statistics, whether to store statistics about local search applications
     */
    template <bool handle_partial_solutions = false, bool log_statistics = false>
    class RandomizedVariableNeighborhoodDescent
        : EnableStatistics<log_statistics, handle_partial_solutions>,
          private NonCopyable<RandomizedVariableNeighborhoodDescent<handle_partial_solutions, log_statistics>>,
          public VariableNeighborhoodDescentInterface {

        const cobra::Instance& instance;
        MoveGenerators& moves;
        std::mt19937& rand_engine;

        // todo CONSTEXPR
        std::unordered_map<Operator, std::function<AbstractOperator*()>> OperatorInitTable;

        std::vector<AbstractOperator*> fixed_list;
        std::vector<AbstractOperator*> sortable_list;

    public:
        RandomizedVariableNeighborhoodDescent(const cobra::Instance& instance_, MoveGenerators& moves_, std::initializer_list<Operator> operators,
                                              std::mt19937& rand_engine_, float tolerance = 0.01f)
            : instance(instance_), moves(moves_), rand_engine(rand_engine_) {

            OperatorInitTable[E10] =    [this, tolerance]() { return new CommonOperator<OneZeroExchange, handle_partial_solutions,  log_statistics>(instance, moves, tolerance); };
            OperatorInitTable[E11] =    [this, tolerance]() { return new CommonOperator<OneOneExchange, handle_partial_solutions,  log_statistics>(instance, moves, tolerance); };
            OperatorInitTable[E20] =    [this, tolerance]() { return new CommonOperator<TwoZeroExchange, handle_partial_solutions,  log_statistics>(instance, moves, tolerance); };
            OperatorInitTable[E21] =    [this, tolerance]() { return new CommonOperator<TwoOneExchange, handle_partial_solutions,  log_statistics>(instance, moves, tolerance); };
            OperatorInitTable[E22] =    [this, tolerance]() { return new CommonOperator<TwoTwoExchange, handle_partial_solutions,  log_statistics>(instance, moves, tolerance); };
            OperatorInitTable[E30] =    [this, tolerance]() { return new CommonOperator<ThreeZeroExchange, handle_partial_solutions,  log_statistics>(instance, moves, tolerance); };
            OperatorInitTable[E31] =    [this, tolerance]() { return new CommonOperator<ThreeOneExchange, handle_partial_solutions,  log_statistics>(instance, moves, tolerance); };
            OperatorInitTable[E32] =    [this, tolerance]() { return new CommonOperator<ThreeTwoExchange, handle_partial_solutions,  log_statistics>(instance, moves, tolerance); };
            OperatorInitTable[E33] =    [this, tolerance]() { return new CommonOperator<ThreeThreeExchange, handle_partial_solutions,  log_statistics>(instance, moves, tolerance); };
            OperatorInitTable[SPLIT] =  [this, tolerance]() { return new CommonOperator<SplitExchange, handle_partial_solutions,  log_statistics>(instance, moves, tolerance); };
            OperatorInitTable[TAILS] =  [this, tolerance]() { return new CommonOperator<TailsExchange, handle_partial_solutions,  log_statistics>(instance, moves, tolerance); };
            OperatorInitTable[TWOPT] =  [this, tolerance]() { return new CommonOperator<TwoOptExchange, handle_partial_solutions,  log_statistics>(instance, moves, tolerance); };
            OperatorInitTable[EJCH] =   [this, tolerance]() { return new CommonOperator<EjectionChain, handle_partial_solutions,  log_statistics, false /*dummy*/,25>(instance, moves, tolerance); };
            OperatorInitTable[RE20] =   [this, tolerance]() { return new CommonOperator<RevTwoZeroExchange, handle_partial_solutions,  log_statistics>(instance, moves, tolerance); };
            OperatorInitTable[RE21] =   [this, tolerance]() { return new CommonOperator<RevTwoOneExchange, handle_partial_solutions,  log_statistics>(instance, moves, tolerance); };
            OperatorInitTable[RE22B] =  [this, tolerance]() { return new CommonOperator<RevTwoTwoExchange, handle_partial_solutions,  log_statistics, true>(instance, moves, tolerance); };
            OperatorInitTable[RE22S] =  [this, tolerance]() { return new CommonOperator<RevTwoTwoExchange, handle_partial_solutions,  log_statistics, false>(instance, moves, tolerance); };
            OperatorInitTable[RE30] =   [this, tolerance]() { return new CommonOperator<RevThreeZeroExchange, handle_partial_solutions,  log_statistics>(instance, moves, tolerance); };
            OperatorInitTable[RE31] =   [this, tolerance]() { return new CommonOperator<RevThreeOneExchange, handle_partial_solutions,  log_statistics>(instance, moves, tolerance); };
            OperatorInitTable[RE32B] =  [this, tolerance]() { return new CommonOperator<RevThreeTwoExchange, handle_partial_solutions,  log_statistics, true>(instance, moves, tolerance); };
            OperatorInitTable[RE32S] =  [this, tolerance]() { return new CommonOperator<RevThreeTwoExchange, handle_partial_solutions,  log_statistics, false>(instance, moves, tolerance); };
            OperatorInitTable[RE33B] =  [this, tolerance]() { return new CommonOperator<RevThreeThreeExchange, handle_partial_solutions,  log_statistics, true>(instance, moves, tolerance); };
            OperatorInitTable[RE33S] =  [this, tolerance]() { return new CommonOperator<RevThreeThreeExchange, handle_partial_solutions,  log_statistics, false>(instance, moves, tolerance); };
            OperatorInitTable[TLCH] =   [this, tolerance]() { return new CommonOperator<TailChain, handle_partial_solutions, log_statistics, false /*dummy*/, 50>(instance, moves, tolerance); };
            OperatorInitTable[STCH] =   [this, tolerance]() { return new CommonOperator<SplitChain, handle_partial_solutions, log_statistics, false /*dummy*/, 50>(instance, moves, tolerance); };

            for (auto op : operators) {
                auto ptr = OperatorInitTable[op]();
                fixed_list.push_back(ptr);
                sortable_list.push_back(ptr);
                if constexpr (log_statistics) {
                    this->operator_ptr_to_id[ptr] = op;
                }
            }
        }

        virtual ~RandomizedVariableNeighborhoodDescent() {
            for (auto move : fixed_list) {
                delete move;
            }
        }

        void apply(cobra::Solution& solution) override {

            std::shuffle(sortable_list.begin(), sortable_list.end(), rand_engine);

            auto end = 0u;

            if constexpr (log_statistics) {
                this->stats.begin(solution);
            }

            auto curr = end;

            do {

                auto ptr = sortable_list[curr];

                const auto improved = ptr->apply_rough_best_improvement(solution);

                if constexpr (log_statistics) {
                    this->stats.append(this->operator_ptr_to_id[sortable_list[curr]], solution);
                }

                if (improved) {
                    end = curr;
                }

                curr = (curr + 1) % sortable_list.size();

            } while (curr != end);

            if constexpr (log_statistics) {
                this->stats.end();
            }


            assert(solution.is_feasible());
        }

        template <typename dummy = Operator>
        auto get_id_from_operator_ptr(AbstractOperator* ptr) -> std::enable_if_t<log_statistics, dummy> {
            return this->operator_ptr_to_id[ptr];
        }

        auto get_operators() -> std::vector<AbstractOperator*>& { return fixed_list; }
    };

    class VariableNeighborhoodDescentComposer {

        const float tolerance;
        std::vector<VariableNeighborhoodDescentInterface*> tiers;

    public:
        VariableNeighborhoodDescentComposer(float tolerance_) : tolerance(tolerance_){};

        void append(VariableNeighborhoodDescentInterface* vnd) { tiers.push_back(vnd); }

        void sequential_apply(Solution& solution) {

        __again__:
            for (auto n = 0u; n < tiers.size(); n++) {
                const auto curr_cost = solution.get_cost();
                tiers[n]->apply(solution);
                if (n > 0 && solution.get_cost() + tolerance < curr_cost) {
                    goto __again__;
                }
            }
        }

        std::vector<VariableNeighborhoodDescentInterface*>& get_tiers() { return tiers; }
    };


}  // namespace cobra


#endif  // COBRA_LOCALSEARCH_HPP