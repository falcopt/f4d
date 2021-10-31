#ifndef _F4D_ABSTRACTOPERATOR_HPP_
#define _F4D_ABSTRACTOPERATOR_HPP_

#include "../Instance.hpp"
#include "../MoveGenerators.hpp"
#include "../Solution.hpp"
#include "../Welford.hpp"

#define UPDATE_BITS_FIRST (0)
#define UPDATE_BITS_SECOND (1)

namespace cobra {

    template <bool log_statistics>
    struct AbstractOperatorStatistics { };
    template <>
    struct AbstractOperatorStatistics<true> {

        // temp variables
        unsigned long long int checked_moves = 0;
        unsigned long long int applied_moves = 0;
        unsigned long long int checked_before_application = 0;
        float old_solution_cost = 0.0f;

        // statistics
        Welford affected_vertices_size;
        Welford num_checked_moves;
        Welford num_applied_moves;
        Welford num_checks_before_application;

        unsigned long long int successful_descent_applications = 0;
        long double total_improvement = 0.0;
    };

    class AbstractOperator : public NonCopyable<AbstractOperator> {
    public:
        AbstractOperator(const Instance& instance_, MoveGenerators& moves_, float tolerance_)
            : instance(instance_), moves(moves_), heap(moves.get_heap()), tolerance(tolerance_), update_bits(moves_.get_update_bits()) { }

        virtual bool apply_rough_best_improvement(Solution& solution) = 0;
        virtual bool apply_best_improvement(Solution& solution) = 0;
        virtual ~AbstractOperator() = default;

    protected:
        virtual void pre_processing(Solution& solution) = 0;
        virtual float compute_cost(const Solution& solution, const MoveGenerator& move) = 0;

        virtual bool is_feasible(const Solution& solution, const MoveGenerator& move) = 0;
        virtual void execute(Solution& solution, const MoveGenerator& move, VertexSet& storage) = 0;
        virtual void post_processing(Solution& solution) = 0;
        virtual std::string get_additional_statistics() = 0;

    protected:
        const Instance& instance;
        MoveGenerators& moves;
        MoveGeneratorsHeap& heap;
        const float tolerance;
        Flat2DVector<bool>& update_bits;
    };

    template <template <auto...> class Tmpl, bool handle_partial_solutions = false, bool log_statistics = false, auto... Args>
    class CommonOperator : public Tmpl<handle_partial_solutions, Args...>, public AbstractOperatorStatistics<log_statistics> {

        using T = Tmpl<handle_partial_solutions, Args...>;

    private:
        TimestampGenerator& timegen;
        VertexSet affected_vertices;


    private:
        inline void symmetric_init(const Solution& solution) {
            const auto currenttime = timegen.get() + 1;
            auto& vtimestamp = T::moves.get_vertex_timestamp();

            auto depot = false;
            for (auto i = solution.get_cache_begin(); i != solution.get_cache_end(); i = solution.get_cache_next(i)) {

                if constexpr (handle_partial_solutions) {
                    if (!solution.is_vertex_in_solution(i)) {
                        continue;
                    }
                }

                // handle the depot as last vertex so that we have better chances to compact move generators and re-use cached info
                if (i == T::instance.get_depot()) {
                    depot = true;
                    continue;
                }

                const auto icache = T::prepare_cache12(solution, i);

                // iterate over 1st indices so that we can exploit the cache
                for (auto move_i1st_index : T::moves.get_move_generator_indices_involving_1st(i)) {

                    const auto& move_i1st = T::moves.get(move_i1st_index);
                    const auto j = move_i1st.get_second_vertex();

                    if constexpr (handle_partial_solutions) {
                        if (!solution.is_vertex_in_solution(j)) {
                            continue;
                        }
                    }

                    if (vtimestamp[j] == currenttime) {  // we already processed j => we already inserted the move gen
                        continue;
                    }

                    // however, consider the base index so that `i` and `j` share the same set of move generators
                    const auto move_index = MoveGenerators::get_base_move_generator_index(move_i1st_index);
                    auto& move = T::moves.get(move_index);

                    const auto jcache = j == T::instance.get_depot() ? T::prepare_cache12(solution, j, i) : T::prepare_cache12(solution, j);

                    const auto delta = T::compute_cost(icache, jcache);
                    move.set_delta(delta);
                    move.set_heap_index(MoveGeneratorsHeap::unheaped);
                    if (move.get_delta() < -T::tolerance) {
                        T::heap.insert(&move);
                    }
                }

                vtimestamp[i] = currenttime;
            }

            if (depot) {

                const auto i = T::instance.get_depot();

                for (auto move_i1st_index : T::moves.get_move_generator_indices_involving_1st(i)) {

                    const auto& move_i1st = T::moves.get(move_i1st_index);
                    const auto j = move_i1st.get_second_vertex();

                    if constexpr (handle_partial_solutions) {
                        if (!solution.is_vertex_in_solution(j)) {
                            continue;
                        }
                    }

                    if (vtimestamp[j] == currenttime) {  // we already processed j => we already inserted the move gen
                        continue;
                    }

                    // however, consider the base index so that `i` and `j` share the same set of move generators
                    const auto move_index = MoveGenerators::get_base_move_generator_index(move_i1st_index);
                    auto& move = T::moves.get(move_index);

                    const auto icache = T::prepare_cache12(solution, i, j);
                    const auto jcache = T::prepare_cache12(solution, j);  // j cannot be the depot

                    const auto delta = T::compute_cost(icache, jcache);
                    move.set_delta(delta);
                    move.set_heap_index(MoveGeneratorsHeap::unheaped);
                    if (move.get_delta() < -T::tolerance) {
                        T::heap.insert(&move);
                    }
                }

                vtimestamp[i] = currenttime;  // todo probably useless
            }

            timegen.increment();
        }

        inline void asymmetric_init(const Solution& solution) {

            const auto currenttime = timegen.get() + 1;

            auto& vtimestamp = T::moves.get_vertex_timestamp();

            auto depot = false;
            for (auto i = solution.get_cache_begin(); i != solution.get_cache_end(); i = solution.get_cache_next(i)) {

                if constexpr (handle_partial_solutions) {
                    if (!solution.is_vertex_in_solution(i)) {
                        continue;
                    }
                }

                // handle the depot as last vertex so that we have better chances to compact move generators and re-use cached info
                if (i == T::instance.get_depot()) {
                    depot = true;
                    continue;
                }

                const auto icache = T::prepare_cache12(solution, i);

                for (const auto move_index : T::moves.get_move_generator_indices_involving_1st(i)) {

                    auto& move = T::moves.get(move_index);

                    const auto j = move.get_second_vertex();

                    if constexpr (handle_partial_solutions) {
                        if (!solution.is_vertex_in_solution(j)) {
                            continue;
                        }
                    }

                    if (vtimestamp[j] == currenttime) {
                        continue;
                    }

                    const auto jcache = j == T::instance.get_depot() ? T::prepare_cache12(solution, j, i) : T::prepare_cache12(solution, j);

                    const auto [delta1, delta2] = T::compute_cost_pair(icache, jcache);

                    move.set_delta(delta1);
                    move.set_heap_index(MoveGeneratorsHeap::unheaped);
                    if (move.get_delta() < -T::tolerance) {
                        T::heap.insert(&move);
                    }

                    const auto twin_move_index = MoveGenerators::get_twin_move_generator_index(move_index);
                    auto& twin_move = T::moves.get(twin_move_index);
                    twin_move.set_delta(delta2);
                    twin_move.set_heap_index(MoveGeneratorsHeap::unheaped);
                    if (twin_move.get_delta() < -T::tolerance) {
                        T::heap.insert(&twin_move);
                    }
                }

                vtimestamp[i] = currenttime;
            }

            if (depot) {

                const auto i = T::instance.get_depot();

                for (const auto move_index : T::moves.get_move_generator_indices_involving_1st(i)) {

                    auto& move = T::moves.get(move_index);

                    const auto j = move.get_second_vertex();

                    if constexpr (handle_partial_solutions) {
                        if (!solution.is_vertex_in_solution(j)) {
                            continue;
                        }
                    }

                    if (vtimestamp[j] == currenttime) {
                        continue;
                    }

                    const auto icache = T::prepare_cache12(solution, i, j);
                    const auto jcache = T::prepare_cache12(solution, j);  // j cannot be the depot

                    const auto [delta1, delta2] = T::compute_cost_pair(icache, jcache);

                    move.set_delta(delta1);
                    move.set_heap_index(MoveGeneratorsHeap::unheaped);
                    if (move.get_delta() < -T::tolerance) {
                        T::heap.insert(&move);
                    }

                    const auto twin_move_index = MoveGenerators::get_twin_move_generator_index(move_index);
                    auto& twin_move = T::moves.get(twin_move_index);
                    twin_move.set_delta(delta2);
                    twin_move.set_heap_index(MoveGeneratorsHeap::unheaped);
                    if (twin_move.get_delta() < -T::tolerance) {
                        T::heap.insert(&twin_move);
                    }
                }

                vtimestamp[i] = currenttime;  // todo probably useless
            }

            timegen.increment();
        }

        inline void initialize_descriptors(const Solution& solution) {

            if constexpr (T::is_symmetric) {
                symmetric_init(solution);
            } else {
                asymmetric_init(solution);
            }
        }

        inline void asymmetric_update(const Solution& solution) {

            const auto currenttime = timegen.get() + 1;

            auto& vtimestamp = T::moves.get_vertex_timestamp();

            const auto heap_insert = [this](MoveGenerator& move, float delta) {
                if (delta > -T::tolerance) {

                    if (move.get_heap_index() != MoveGeneratorsHeap::unheaped) {
                        T::heap.remove(move.get_heap_index());
                    }

                    move.set_delta(delta);
                } else {

                    if (move.get_heap_index() == MoveGeneratorsHeap::unheaped) {
                        move.set_delta(delta);
                        T::heap.insert(&move);
                    } else {
                        T::heap.change_value(move.get_heap_index(), delta);
                    }
                }
            };

            auto depot = false;
            for (auto i : affected_vertices.get_vertices()) {

                if constexpr (handle_partial_solutions) {
                    if (!solution.is_vertex_in_solution(i)) {
                        continue;
                    }
                }

                // handle the depot as last vertex so that we have better chances to compact move generators and re-use cached info
                if (i == T::instance.get_depot()) {
                    depot = true;
                    continue;
                }

                const auto iupij = T::update_bits.at(i, UPDATE_BITS_FIRST);
                const auto iupji = T::update_bits.at(i, UPDATE_BITS_SECOND);

                if (iupij && iupji) {  // update both (i, j) and (j, i)

                    const auto icache = T::prepare_cache12(solution, i);

                    for (const auto move_index : T::moves.get_move_generator_indices_involving_1st(i)) {

                        auto& move = T::moves.get(move_index);
                        const auto j = move.get_second_vertex();

                        if constexpr (handle_partial_solutions) {
                            if (!solution.is_vertex_in_solution(j)) {
                                continue;
                            }
                        }

                        if (vtimestamp[j] == currenttime) {  // move gen (j, i) and (i, j) may have been already updated

                            const auto jupji = T::update_bits.at(j, UPDATE_BITS_FIRST);
                            const auto jupij = T::update_bits.at(j, UPDATE_BITS_SECOND);

                            if (jupji && jupij) {
                            }                  // (i, j) and (j, i) were already updated
                            else if (jupji) {  // (j, i) already updated => update (i, j) only

                                const auto jcache = j == T::instance.get_depot() ? T::prepare_cache2(solution, j, i) : T::prepare_cache2(solution, j);
                                const auto delta = T::compute_cost(icache, jcache);

                                heap_insert(move, delta);

                            } else if (jupij) {  // (i, j) already updated => update (j, i) only

                                const auto jcache = j == T::instance.get_depot() ? T::prepare_cache1(solution, j, i) : T::prepare_cache1(solution, j);
                                const auto twin_delta = T::compute_cost(jcache, icache);

                                const auto twin_move_index = MoveGenerators::get_twin_move_generator_index(move_index);
                                auto& twin_move = T::moves.get(twin_move_index);

                                heap_insert(twin_move, twin_delta);
                            }

                        } else {  // j was not updated before

                            const auto jcache = j == T::instance.get_depot() ? T::prepare_cache12(solution, j, i) : T::prepare_cache12(solution, j);
                            const auto [delta1, delta2] = T::compute_cost_pair(icache, jcache);

                            const auto twin_move_index = MoveGenerators::get_twin_move_generator_index(move_index);
                            auto& twin_move = T::moves.get(twin_move_index);

                            heap_insert(move, delta1);
                            heap_insert(twin_move, delta2);
                        }
                    }

                } else if (iupij) {  // update only (i, j)

                    const auto icache = T::prepare_cache1(solution, i);

                    for (const auto move_index : T::moves.get_move_generator_indices_involving_1st(i)) {

                        auto& move = T::moves.get(move_index);
                        const auto j = move.get_second_vertex();

                        if constexpr (handle_partial_solutions) {
                            if (!solution.is_vertex_in_solution(j)) {
                                continue;
                            }
                        }

                        if (vtimestamp[j] != currenttime ||  // j has not already been processed
                            (vtimestamp[j] == currenttime &&
                             !T::update_bits.at(j, UPDATE_BITS_SECOND)))  // j was processed but (i, j) was not updated because not required
                        {

                            const auto jcache = j == T::instance.get_depot() ? T::prepare_cache2(solution, j, i) : T::prepare_cache2(solution, j);
                            const auto delta = T::compute_cost(icache, jcache);

                            heap_insert(move, delta);
                        }
                    }

                } else if (iupji) {  // update only (j, i)

                    const auto icache = T::prepare_cache2(solution, i);

                    for (const auto move_index : T::moves.get_move_generator_indices_involving_2nd(i)) {

                        auto& move = T::moves.get(move_index);
                        const auto j = move.get_first_vertex();

                        if constexpr (handle_partial_solutions) {
                            if (!solution.is_vertex_in_solution(j)) {
                                continue;
                            }
                        }

                        if (vtimestamp[j] != currenttime ||  // j has not already been processed
                            (vtimestamp[j] == currenttime &&
                             !T::update_bits.at(j, UPDATE_BITS_FIRST)))  // j was processed but (j, i) was not updated because not required
                        {

                            const auto jcache = j == T::instance.get_depot() ? T::prepare_cache1(solution, j, i) : T::prepare_cache1(solution, j);
                            const auto delta = T::compute_cost(jcache, icache);

                            heap_insert(move, delta);
                        }
                    }
                }

                vtimestamp[i] = currenttime;
            }

            if (depot) {

                const auto i = T::instance.get_depot();

                const auto iupij = T::update_bits.at(i, UPDATE_BITS_FIRST);
                const auto iupji = T::update_bits.at(i, UPDATE_BITS_SECOND);

                if (iupij && iupji) {  // update both (i, j) and (j, i)

                    for (const auto move_index : T::moves.get_move_generator_indices_involving_1st(i)) {

                        auto& move = T::moves.get(move_index);

                        const auto j = move.get_second_vertex();

                        if constexpr (handle_partial_solutions) {
                            if (!solution.is_vertex_in_solution(j)) {
                                continue;
                            }
                        }

                        if (vtimestamp[j] == currenttime) {  // move gen (j, i) and (i, j) may have been already updated

                            const auto jupji = T::update_bits.at(j, UPDATE_BITS_FIRST);
                            const auto jupij = T::update_bits.at(j, UPDATE_BITS_SECOND);

                            if (jupji && jupij) {
                            }                  // (i, j) and (j, i) were already updated
                            else if (jupji) {  // (j, i) already updated

                                const auto icache = T::prepare_cache1(solution, i, j);
                                const auto jcache = T::prepare_cache2(solution, j);
                                const auto delta = T::compute_cost(icache, jcache);

                                heap_insert(move, delta);

                            } else if (jupij) {  // (i, j) already updated

                                const auto twin_move_index = MoveGenerators::get_twin_move_generator_index(move_index);
                                auto& twin_move = T::moves.get(twin_move_index);

                                const auto icache = T::prepare_cache2(solution, i, j);
                                const auto jcache = T::prepare_cache1(solution, j);
                                const auto twin_delta = T::compute_cost(jcache, icache);

                                heap_insert(twin_move, twin_delta);
                            }

                        } else {  // j was not updated before

                            const auto icache = T::prepare_cache12(solution, i, j);
                            const auto jcache = T::prepare_cache12(solution, j);
                            const auto [delta1, delta2] = T::compute_cost_pair(icache, jcache);

                            heap_insert(move, delta1);

                            const auto twin_move_index = MoveGenerators::get_twin_move_generator_index(move_index);
                            auto& twin_move = T::moves.get(twin_move_index);
                            heap_insert(twin_move, delta2);
                        }
                    }

                } else if (iupij) {  // update only (i, j)


                    for (const auto move_index : T::moves.get_move_generator_indices_involving_1st(i)) {

                        auto& move = T::moves.get(move_index);
                        const auto j = move.get_second_vertex();

                        if constexpr (handle_partial_solutions) {
                            if (!solution.is_vertex_in_solution(j)) {
                                continue;
                            }
                        }

                        if (vtimestamp[j] != currenttime ||  // j has not already been processed
                            (vtimestamp[j] == currenttime &&
                             !T::update_bits.at(j, UPDATE_BITS_SECOND)))  // j was processed but (i, j) was not updated because not required
                        {

                            const auto icache = T::prepare_cache1(solution, i, j);
                            const auto jcache = T::prepare_cache2(solution, j);
                            const auto delta = T::compute_cost(icache, jcache);
                            heap_insert(move, delta);
                        }
                    }

                } else if (iupji) {  // update only (j, i)

                    for (const auto move_index : T::moves.get_move_generator_indices_involving_2nd(i)) {

                        auto& move = T::moves.get(move_index);
                        const auto j = move.get_first_vertex();

                        if constexpr (handle_partial_solutions) {
                            if (!solution.is_vertex_in_solution(j)) {
                                continue;
                            }
                        }

                        if (vtimestamp[j] != currenttime ||  // j has not already been processed
                            (vtimestamp[j] == currenttime &&
                             !T::update_bits.at(j, UPDATE_BITS_FIRST)))  // j was processed but (j, i) was not updated because not required
                        {

                            const auto icache = T::prepare_cache2(solution, i, j);
                            const auto jcache = T::prepare_cache1(solution, j);
                            const auto delta = T::compute_cost(jcache, icache);
                            heap_insert(move, delta);
                        }
                    }
                }

                vtimestamp[i] = currenttime;  // todo probably useless
            }

            timegen.increment();
        }

        inline void symmetric_update(const Solution& solution) {

            const auto currenttime = timegen.get() + 1;

            auto& vtimestamp = T::moves.get_vertex_timestamp();

            const auto heap_insert = [this](MoveGenerator& move, float delta) {
                if (delta > -T::tolerance) {

                    if (move.get_heap_index() != MoveGeneratorsHeap::unheaped) {
                        T::heap.remove(move.get_heap_index());
                    }

                    move.set_delta(delta);
                } else {

                    if (move.get_heap_index() == MoveGeneratorsHeap::unheaped) {
                        move.set_delta(delta);
                        T::heap.insert(&move);
                    } else {
                        T::heap.change_value(move.get_heap_index(), delta);
                    }
                }
            };

            auto depot = false;
            for (auto i : affected_vertices.get_vertices()) {

                if constexpr (handle_partial_solutions) {
                    if (!solution.is_vertex_in_solution(i)) {
                        continue;
                    }
                }

                // handle the depot as last vertex so that we have better chances to compact move generators and re-use cached info
                if (i == T::instance.get_depot()) {
                    depot = true;
                    continue;
                }

                auto icache = T::prepare_cache12(solution, i);

                // iterate over 1st indices so that we can exploit the cache
                for (auto move_i1st_index : T::moves.get_move_generator_indices_involving_1st(i)) {

                    const auto& move_i1st = T::moves.get(move_i1st_index);
                    const auto j = move_i1st.get_second_vertex();

                    if constexpr (handle_partial_solutions) {
                        if (!solution.is_vertex_in_solution(j)) {
                            continue;
                        }
                    }

                    if (vtimestamp[j] == currenttime) {  // we already processed j => we already considered the move gen
                        continue;
                    }

                    // however, consider the base index so that `i` and `j` share the same set of move generators
                    const auto move_index = MoveGenerators::get_base_move_generator_index(move_i1st_index);
                    auto& move = T::moves.get(move_index);

                    const auto jcache = j == T::instance.get_depot() ? T::prepare_cache12(solution, j, i) : T::prepare_cache12(solution, j);

                    const auto delta = T::compute_cost(icache, jcache);
                    heap_insert(move, delta);
                }
            }

            if (depot) {

                const auto i = T::instance.get_depot();

                for (auto move_i1st_index : T::moves.get_move_generator_indices_involving_1st(i)) {

                    const auto& move_i1st = T::moves.get(move_i1st_index);
                    const auto j = move_i1st.get_second_vertex();

                    if constexpr (handle_partial_solutions) {
                        if (!solution.is_vertex_in_solution(j)) {
                            continue;
                        }
                    }

                    if (vtimestamp[j] == currenttime) {  // we already processed j => we already inserted the move gen
                        continue;
                    }

                    // however, consider the base index so that `i` and `j` share the same set of move generators
                    const auto move_index = MoveGenerators::get_base_move_generator_index(move_i1st_index);
                    auto& move = T::moves.get(move_index);

                    const auto icache = T::prepare_cache12(solution, i, j);
                    const auto jcache = T::prepare_cache12(solution, j);  // j cannot be the depot

                    const auto delta = T::compute_cost(icache, jcache);
                    heap_insert(move, delta);
                }

                vtimestamp[i] = currenttime;  // todo probably useless
            }

            timegen.increment();
        }

        inline void descriptors_update(const Solution& solution) {

            if constexpr (T::is_symmetric) {
                symmetric_update(solution);
            } else {
                asymmetric_update(solution);
            }
        }

    public:
        CommonOperator(const Instance& instance_, MoveGenerators& moves_, float tolerance_)
            : T(instance_, moves_, tolerance_), timegen(moves_.get_timestamp_generator()), affected_vertices(instance_.get_vertices_num()) { }

        bool apply_rough_best_improvement(Solution& solution) {

            T::heap.reset();

            T::pre_processing(solution);

            initialize_descriptors(solution);

            if constexpr (log_statistics) {
                this->checked_moves = 0;
                this->applied_moves = 0;
                this->checked_before_application = 0;
                this->old_solution_cost = solution.get_cost();
            }

            auto improved = false;

            auto index = 0;

            while (index < T::heap.size()) {
                auto& move = *T::heap.spy(index++);

                if constexpr (log_statistics) {
                    this->checked_moves++;
                    this->checked_before_application++;
                }

                if constexpr (handle_partial_solutions) {
                    if (!solution.is_vertex_in_solution(move.get_first_vertex()) || !solution.is_vertex_in_solution(move.get_second_vertex())) {
                        continue;
                    }
                }

                if (!T::is_feasible(solution, move)) {
                    continue;
                }

#ifndef NDEBUG
                auto old_cost = solution.get_cost();
#endif

                T::execute(solution, move, affected_vertices);

                if constexpr (log_statistics) {
                    this->applied_moves++;
                    this->affected_vertices_size.update(affected_vertices.size());
                    this->num_checks_before_application.update(this->checked_before_application);
                    this->checked_before_application = 0;
                }

                // assert(std::fabs(solution.get_cost() - old_cost - move.get_delta()) < T::tolerance); wrong for ejection chain!
                assert(old_cost > solution.get_cost());

                assert(solution.is_feasible());

                improved = true;
                index = 0;

                descriptors_update(solution);

                affected_vertices.clear();
            }

            T::post_processing(solution);

            if constexpr (log_statistics) {
                this->num_applied_moves.update(this->applied_moves);
                this->num_checked_moves.update(this->checked_moves);
                if (improved) {
                    this->successful_descent_applications++;
                    this->total_improvement += (this->old_solution_cost - solution.get_cost());
                }
            }

            return improved;
        }

        bool apply_best_improvement(Solution& solution) {

            T::heap.reset();

            T::pre_processing(solution);

            initialize_descriptors(solution);

            auto improved = false;

            auto infeasible_smds = std::vector<MoveGenerator*>();

            while (!T::heap.is_empty()) {

                auto& move = *T::heap.get();

                if constexpr (handle_partial_solutions) {
                    if (!solution.is_vertex_in_solution(move.get_first_vertex()) || !solution.is_vertex_in_solution(move.get_second_vertex())) {
                        continue;
                    }
                }

                if (!T::is_feasible(solution, move)) {
                    infeasible_smds.push_back(&move);
                    continue;
                }

#ifndef NDEBUG
                auto old_cost = solution.get_cost();
#endif

                T::execute(solution, move, affected_vertices);

                assert(old_cost > solution.get_cost());

                assert(solution.is_feasible());

                improved = true;

                for (auto smd : infeasible_smds) {
                    T::heap.insert(smd);
                }
                infeasible_smds.clear();

                descriptors_update(solution);

                affected_vertices.clear();
            }

            return improved;
        }

        template <typename dummy = std::string>
        auto get_statistics() -> std::enable_if_t<log_statistics, dummy> {
            std::string text;

            text.append("MeanNumAffectedVertices \t " + std::to_string(this->affected_vertices_size.get_mean()) + "\n");
            text.append("StdDevNumAffectedVertices \t " + std::to_string(this->affected_vertices_size.get_standard_deviation()) + "\n");

            text.append("MeanNumCheckedMoves \t " + std::to_string(this->num_checked_moves.get_mean()) + "\n");
            text.append("StdDevNumCheckedMoves \t " + std::to_string(this->num_checked_moves.get_standard_deviation()) + "\n");

            text.append("MeanNumAppliedMoves \t " + std::to_string(this->num_applied_moves.get_mean()) + "\n");
            text.append("StdDevNumAppliedMoves \t " + std::to_string(this->num_applied_moves.get_standard_deviation()) + "\n");

            text.append("MeanNumCheckedMovesBeforeApplication \t " + std::to_string(this->num_checks_before_application.get_mean()) + "\n");
            text.append("StdDevNumCheckedMovesBeforeApplication \t " + std::to_string(this->num_checks_before_application.get_standard_deviation()) + "\n");

            text.append("TotalOperatorImprovement \t " + std::to_string(this->total_improvement) + "\n");
            text.append("SuccessfulDescentApplications \t " + std::to_string(this->successful_descent_applications) + "\n");

            text.append(this->get_additional_statistics());

            return text;
        }
    };


}  // namespace cobra
#endif