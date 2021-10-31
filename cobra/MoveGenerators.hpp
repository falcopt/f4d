#ifndef _F4D_MOVEGENERATORS_HPP_
#define _F4D_MOVEGENERATORS_HPP_

#include <functional>
#include <random>
#include <set>
#include <unordered_set>

#include "BinaryHeap.hpp"
#include "Instance.hpp"
#include "VectorView.hpp"
#include "VertexSet.hpp"

namespace cobra {

    /**
     * Simple counter
     */
    class TimestampGenerator : private NonCopyable<TimestampGenerator> {
        unsigned long value = 0;

    public:
        TimestampGenerator() = default;
        inline unsigned long get() {
            return value;
        }
        inline void increment() {
            ++value;
        }
    };

    /**
     * Single move generator / static move descriptor
     */
    class MoveGenerator : private NonCopyable<MoveGenerator> {
        int first_vertex;
        int second_vertex;
        float delta = 0.0f;
        int heap_index = -1;

    public:
        MoveGenerator(int i, int j) : first_vertex(i), second_vertex(j) { }

        inline auto get_first_vertex() const {
            return first_vertex;
        }
        inline auto get_second_vertex() const {
            return second_vertex;
        }

        inline auto get_delta() const {
            return delta;
        }
        inline auto set_delta(float value) {
            delta = value;
        }

        inline auto get_heap_index() const {
            return heap_index;
        }
        inline auto set_heap_index(int index) {
            heap_index = index;
        }
    };


    struct MGCompare {
        auto operator()(MoveGenerator* mg1, MoveGenerator* mg2) {
            assert(mg1 && mg2);
            return mg1->get_delta() - mg2->get_delta();
        }
    };

    struct MGGetIdx {
        auto operator()(MoveGenerator* mg1) {
            assert(mg1);
            return mg1->get_heap_index();
        }
    };

    struct MGSetIdx {
        void operator()(MoveGenerator* mg1, int idx) {
            assert(mg1);
            mg1->set_heap_index(idx);
        }
    };

    struct MGUpdate {
        auto operator()(MoveGenerator* mg1, float delta) {
            assert(mg1);
            const auto res = mg1->get_delta() - delta;
            mg1->set_delta(delta);
            return res;
        }
    };

    class MoveGeneratorsHeap : private NonCopyable<MoveGeneratorsHeap>, private BinaryHeap<MoveGenerator*, MGCompare, MGGetIdx, MGSetIdx, MGUpdate, -1> {

        typedef BinaryHeap<MoveGenerator*, MGCompare, MGGetIdx, MGSetIdx, MGUpdate> BHeap;

        void dump() override {
            for (auto n = 0; n < size(); n++) {
                const auto& move = spy(n);
                std::cout << "[" << n << "] (" << move->get_first_vertex() << ", " << move->get_second_vertex() << ") delta = " << move->get_delta()
                          << " heap index = " << move->get_heap_index() << "\n";
            }
        }

    public:
        MoveGeneratorsHeap() = default;
        MoveGeneratorsHeap(const MoveGeneratorsHeap& other) = default;


        void reset() {
            BHeap::reset();
        }
        bool is_empty() const {
            return BHeap::empty();
        }
        void insert(MoveGenerator* mg) {
            BHeap::insert(mg);
        }
        MoveGenerator* get() {
            return BHeap::get();
        }
        void remove(int heap_index) {
            BHeap::remove(heap_index);
        }
        void change_value(int heap_index, float value) {
            BHeap::update(heap_index, value);
        };
        int size() const {
            return BHeap::size();
        };
        MoveGenerator* spy(int heap_index) {
            return BHeap::spy(heap_index);
        }

        static const int unheaped = -1;
    };


    class MoveGenerators;
    /**
     * Generic view over a set of move generators
     */
    class AbstractMoveGeneratorsView : private NonCopyable<AbstractMoveGeneratorsView> {

        std::function<std::vector<int>(int)> generator;

    protected:
        const cobra::Instance& instance;
        std::vector<std::vector<int>> move_generator_indices_involving;
        std::vector<std::vector<int>> all_move_generator_indices_involving;

        std::unordered_map<int, int> move_map;
        std::vector<std::pair<bool, bool>> active_in;

        void setup_active_trackers() {
            // Each view must keep track of the vertices in which each move generator is active
            for (auto i = instance.get_vertices_begin(); i < instance.get_vertices_end(); i++) {
                for (const auto idx : all_move_generator_indices_involving[i]) {
                    if (move_map.count(idx)) {
                        continue;
                    }
                    move_map.insert({idx, active_in.size()});
                    active_in.emplace_back(false, false);  // at the beginning the move gen is not active in any vertex
                }
            }
        }

        inline auto set_active_in(const MoveGenerator& move, int move_idx, int vertex) {
            assert(move_map.count(move_idx));
            const auto mapped_idx = move_map[move_idx];
            assert(vertex == move.get_first_vertex() || vertex == move.get_second_vertex());
            if (vertex == move.get_first_vertex()) {
                active_in[mapped_idx].first = true;
            } else {
                active_in[mapped_idx].second = true;
            }
        }

        inline auto set_not_active_in(const MoveGenerator& move, int move_idx, int vertex) {
            assert(move_map.count(move_idx));
            const auto mapped_idx = move_map[move_idx];
            assert(vertex == move.get_first_vertex() || vertex == move.get_second_vertex());
            if (vertex == move.get_first_vertex()) {
                active_in[mapped_idx].first = false;
            } else {
                active_in[mapped_idx].second = false;
            }
        }

        inline auto is_active_in(const MoveGenerator& move, int move_idx, int vertex) {
            assert(move_map.count(move_idx));
            const auto mapped_idx = move_map[move_idx];
            assert(vertex == move.get_first_vertex() || vertex == move.get_second_vertex());
            if (vertex == move.get_first_vertex()) {
                return active_in[mapped_idx].first;
            } else {
                return active_in[mapped_idx].second;
            }
        }

        inline auto is_active_in_other(const MoveGenerator& move, int move_idx, int vertex) {
            assert(move_map.count(move_idx));
            const auto mapped_idx = move_map[move_idx];
            assert(vertex == move.get_first_vertex() || vertex == move.get_second_vertex());
            if (vertex == move.get_first_vertex()) {
                return active_in[mapped_idx].second;
            } else {
                return active_in[mapped_idx].first;
            }
        }

    public:
        AbstractMoveGeneratorsView(const cobra::Instance& instance_, std::function<std::vector<int>(int)> generator_)
            : generator(std::move(generator_)),  // `generator_` is no more valid after `std::move`
              instance(instance_) {

            move_generator_indices_involving.resize(instance.get_vertices_num());
            all_move_generator_indices_involving.resize(instance.get_vertices_num());
        }

        /**
         * Returns the generator
         * @return generator
         */
        auto get_generator() -> std::function<std::vector<int>(int)>& {
            return generator;
        }
        auto get_move_generator_indices_involving(int vertex) -> std::vector<int>& {
            return move_generator_indices_involving[vertex];
        }

        /**
         * During `MoveGenerators` construction, each move generator generated by `generator` is assigned a unique index
         * within the `moves` vector in `MoveGenerators`. `__notify_move_indices` is called during this process to
         * notify the view about the indices associated with move generators involving `vertex`
         * @param vertex processing move generators involving `vertex`
         * @param indices list of unique indices for move generators involving `vertex`
         */
        virtual auto __notify_move_indices(int vertex, std::vector<int>& indices, MoveGenerators& moves) -> void = 0;
        /**
         * Notify the view that the `MoveGenerators` object has now completely built the set of move
         * generators. Views can use this callback to perform post-processing actions.
         * @param moves
         */
        virtual auto __notify_build_complete(MoveGenerators& moves) -> void = 0;
        virtual auto set_active_percentage(std::vector<float>& percentage, std::vector<int>& vertices, MoveGenerators& moves,
                                           cobra::VertexSet& vertices_in_updated_moves) -> void = 0;
        auto get_number_of_moves() -> int {

            auto set = std::set<int>();
            for (auto i = instance.get_vertices_begin(); i < instance.get_vertices_end(); i++) {
                for (auto id : all_move_generator_indices_involving[i]) {
                    set.insert(id);
                }
            }
            return 2 * static_cast<int>(set.size());  // here we are storing a single copy between (i, j) and (j,i)
        }

        inline auto is_active(int move_idx) {
            assert(move_map.count(move_idx));
            const auto mapped_idx = move_map[move_idx];
            return active_in[mapped_idx].first || active_in[mapped_idx].second;
        }
    };

    /**
     * Move generators container
     */
    class MoveGenerators : private NonCopyable<MoveGenerators> {

        const cobra::Instance& instance;
        std::vector<MoveGenerator> moves;
        std::vector<AbstractMoveGeneratorsView*>& views;
        MoveGeneratorsHeap heap;
        std::vector<std::vector<int>> move_generator_indices_involving;  // summarize view indices without duplicates
        std::vector<float> prev_percentage;

        TimestampGenerator timegen;
        Flat2DVector<bool> update_bits;

        std::vector<unsigned long> vertex_timestamp;

    public:
        MoveGenerators(const cobra::Instance& instance_, std::vector<AbstractMoveGeneratorsView*>& views_)
            : instance(instance_), views(views_), heap(MoveGeneratorsHeap()), vertex_timestamp(instance_.get_vertices_num(), 0) {

            update_bits.resize(instance.get_vertices_num(), 2);

            struct pair_hash {
                auto operator()(const std::pair<int, int>& p) const -> size_t {
                    const auto prime = 31;
                    auto result = 1;
                    result = prime * result + p.first;
                    result = prime * result + p.second;
                    return std::hash<int>()(result);
                }
            };

            // map move generator -> unique index
            auto unique_moves = std::unordered_map<std::pair<int, int>, int, pair_hash>();

            // identify the set of unique move generators across all views
            for (auto& view : views) {

                const auto& generator = view->get_generator();

                for (auto i = instance.get_vertices_begin(); i < instance.get_vertices_end(); i++) {

                    const auto& endpoints = generator(i);

                    auto move_indices_in_moves = std::vector<int>();
                    move_indices_in_moves.reserve(endpoints.size());

                    for (auto j : endpoints) {

                        auto a = i;
                        auto b = j;
                        if (b < a) {
                            std::swap(a, b);
                        }

                        if (unique_moves.count({a, b})) {

                            const auto move_index = unique_moves[{a, b}];

                            move_indices_in_moves.emplace_back(move_index);

                        } else {

                            const auto move_index = moves.size();
                            move_indices_in_moves.emplace_back(move_index);

                            moves.emplace_back(a, b);
                            moves.emplace_back(b, a);

                            unique_moves[{a, b}] = move_index;
                        }
                    }

                    view->__notify_move_indices(i, move_indices_in_moves, *this);
                }

                view->__notify_build_complete(*this);
            }

            // note that the heap MUST be initialized once all the move generators have been placed into the `moves` vector
            // this->heap = new MoveGeneratorsHeap(*this); //#BH_CHANGE


            this->move_generator_indices_involving.resize(instance.get_vertices_num());

            this->prev_percentage.resize(instance.get_vertices_num(), 0.0f);
        }

        virtual ~MoveGenerators() { }
        /**
         * Returns the move generator indexed `idx`
         * @param index
         * @return move generator indexed `idx`
         */
        inline MoveGenerator& get(int idx) {
            return moves[idx];
        }
        /**
         * Set the percentage of move generators to consider in local search procedures
         * @param percentage float value in [0, 1]
         */
        void set_active_percentage(std::vector<float>& percentage, std::vector<int>& vertices) {

            // Update the active move generators of each view associated with the list of vertices in input

            // First identify all vertices that needs an update
            for (auto n = 0u; n < vertices.size();) {
                const auto vertex = vertices[n];
                if (std::fabs(percentage[vertex] - prev_percentage[vertex]) < 0.01f) {
                    // Remove a vertex if it does not need to be updated
                    std::swap(vertices[n], vertices[vertices.size() - 1]);
                    vertices.pop_back();
                } else {
                    n++;
                }
            }

            if (vertices.empty()) {
                return;
            }

            // For each view, set or unset the move generators involving each 'vertex' according to percentage[vertex]
            // Collect all vertices involving moves that are set or unset
            // auto vertices_in_updated_moves = std::unordered_set<int>();
            auto vertices_in_updated_moves = cobra::VertexSet(instance.get_vertices_num());
            for (auto& view : views) {
                view->set_active_percentage(percentage, vertices, *this, vertices_in_updated_moves);
            }

            // The previous view->set_active_percentage(...) for vertex 'i' may cause a
            // - direct updates: we are adding or removing move generators (i, j) to the list of move generators involving 'i'
            // - indirect updates: we are adding move generators (i, j) to the list of move generators involving 'j' because
            //   of the direct updates (addition) or we are setting as not active move generators (i, j) that were in the list
            //   of move generators involving 'j' (removal). For the latter, when (i, j) was only active due to 'i' and should now
            //   be removed from the list of move generators involving 'j' because no longer active in 'j', we do not have an efficient
            //   way to remove it from there. A possibility would be to store in each view a matrix of move generator positions.
            //   A simpler way is to check and not add move generators that are not active in both vertices.
            //   Note that this latter approach causes the views to always stay in possibly inconsistent states. One can
            //   access the correct list of move generators associated to a given vertex 'i' by accessing the filtered list in
            //   the 'MoveGenerators' class.
            // Here in the following we update those lists.
            // auto unique_move_generator_set = std::unordered_set<int>();
            auto unique_move_generators = std::vector<int>();
            auto unique_endpoints = cobra::VertexSet(instance.get_vertices_num());
            for (const auto vertex : vertices_in_updated_moves.get_vertices()) {
                // unique_move_generator_set.clear();
                unique_move_generators.clear();
                unique_endpoints.clear();
                for (auto& view : views) {
                    for (auto move_idx : view->get_move_generator_indices_involving(vertex)) {
                        if (view->is_active(move_idx)) {

                            if (vertex != moves[move_idx].get_first_vertex()) {
                                move_idx = get_twin_move_generator_index(move_idx);
                            }

                            const auto& move = moves[move_idx];
                            int other_vertex = move.get_second_vertex();

                            if (!unique_endpoints.contains(other_vertex)) {
                                unique_endpoints.insert_without_checking_existance(other_vertex);
                                unique_move_generators.push_back(move_idx);
                            }
                            // unique_move_generator_set.insert(move_idx);
                        }
                    }
                }

                move_generator_indices_involving[vertex] = unique_move_generators;

                // move_generator_indices_involving[vertex].clear();
                // move_generator_indices_involving[vertex].reserve(unique_move_generator_set.size());
                // for(auto move : unique_move_generator_set) { move_generator_indices_involving[vertex].emplace_back(move); }
                // sort indices to obtain the same vector across different compilers/platforms
                // This is unfortunately pretty useless and time consuming.
                // In fact, I forgot that also in LocalSearch.hpp I am iterating over an unordered_set when scanning the vertices
                // affected by a move application...
                // I hate you unordered data structures è.é
                // std::sort(move_generator_indices_involving[vertex].begin(), move_generator_indices_involving[vertex].end());
            }

            // Finally, update the stored percentages
            for (const auto vertex : vertices) {
                prev_percentage[vertex] = percentage[vertex];
            }
        }
        void set_active_percentage(float percentage) {

            auto vertices = std::vector<int>(instance.get_vertices_num());
            auto pi = std::vector<float>(instance.get_vertices_num());

            for (auto i = instance.get_vertices_begin(); i < instance.get_vertices_end(); i++) {
                vertices[i] = i;
                pi[i] = percentage;
            }

            set_active_percentage(pi, vertices);
        }
        /**
         * Returns the heap used for storing move generators during local search applications
         * @return heap
         */
        inline MoveGeneratorsHeap& get_heap() {
            return heap;
        }
        /**
         * Returns the vector of move generators
         * @return vector of move generators
         */
        inline std::vector<MoveGenerator>& get_raw_vector() {
            return moves;
        }

        static inline int get_twin_move_generator_index(int index) {
            return index ^ 1;
        }

        static inline auto get_base_move_generator_index(int index) {
            return index & ~1;
        }

        inline const auto& get_move_generator_indices_involving_1st(int vertex) const {
            return move_generator_indices_involving[vertex];
        }

        inline auto get_move_generator_indices_involving_2nd(int vertex) const {
            const auto& v = move_generator_indices_involving[vertex];
            return VectorView<decltype(v.begin()), twin_functor>(v.begin(), v.end());
        }

        inline auto get_move_generator_indices_involving(int vertex) const {
            const auto& v = move_generator_indices_involving[vertex];
            return VectorView<decltype(v.begin()), base_functor>(v.begin(), v.end());
        }

        inline TimestampGenerator& get_timestamp_generator() {
            return timegen;
        }

        inline Flat2DVector<bool>& get_update_bits() {
            return update_bits;
        }

        inline std::vector<unsigned long>& get_vertex_timestamp() {
            return vertex_timestamp;
        }
    };


    class KNeighborsMoveGeneratorsView : public AbstractMoveGeneratorsView {
        const int max_neighbors_num;          // k max, `set_active_percentage` can be used to vary 0 <= k <= k max
        std::vector<int> curr_neighbors_num;  // current k


        void __notify_move_indices([[maybe_unused]] int vertex, std::vector<int>& indices, MoveGenerators& moves) override {
            for (const auto index : indices) {
                const auto& move = moves.get(index);
                const auto i = move.get_first_vertex();
                const auto j = move.get_second_vertex();
                assert(i == vertex || j == vertex);
                all_move_generator_indices_involving[i].emplace_back(index);
                all_move_generator_indices_involving[j].emplace_back(index);
            }
        }
        void __notify_build_complete(MoveGenerators& moves) override {

            for (auto i = instance.get_customers_begin(); i < instance.get_customers_end(); i++) {

                auto& indices = all_move_generator_indices_involving[i];

                auto set = std::unordered_set<int>();

                set.insert(indices.begin(), indices.end());
                indices.clear();
                indices.insert(indices.begin(), set.begin(), set.end());
                indices.shrink_to_fit();
                std::sort(indices.begin(), indices.end(), [&moves, this](auto a, auto b) {
                    const auto& a_move = moves.get(a);
                    const auto a_cost = this->instance.get_cost(a_move.get_first_vertex(), a_move.get_second_vertex());

                    const auto& b_move = moves.get(b);
                    const auto b_cost = this->instance.get_cost(b_move.get_first_vertex(), b_move.get_second_vertex());

                    return a_cost < b_cost;
                });
            }

            setup_active_trackers();
        }
        void set_active_percentage(std::vector<float>& percentage, std::vector<int>& vertices, MoveGenerators& moves,
                                   cobra::VertexSet& vertices_in_updated_moves) override {

            auto vertices_to_update = std::vector<int>();

            // First, for each vertex identify what are the move generators to set as active or not active
            for (const auto vertex : vertices) {

                const auto new_num = static_cast<int>(std::round(percentage[vertex] * static_cast<float>(max_neighbors_num)));

                assert(new_num <= static_cast<int>(all_move_generator_indices_involving[vertex].size()));

                if (new_num == curr_neighbors_num[vertex]) {
                    continue;
                }

                vertices_to_update.push_back(vertex);

                // Set as active or not active the move generators associated to vertex according to 'new_num'
                // Keep track of the moves involved in the update by storing the vertices in 'vertices_in_updated_moves'
                // This set will be later used by the 'MoveGenerators' class to perform a selective update of the necessary
                // data structures

                if (new_num < curr_neighbors_num[vertex]) {
                    // removal
                    for (auto n = new_num; n < curr_neighbors_num[vertex]; n++) {
                        const auto move_idx = all_move_generator_indices_involving[vertex][n];
                        const auto& move = moves.get(move_idx);
                        this->set_not_active_in(move, move_idx, vertex);
                        vertices_in_updated_moves.insert(move.get_first_vertex());
                        vertices_in_updated_moves.insert(move.get_second_vertex());
                    }
                } else {
                    // addition
                    for (auto n = curr_neighbors_num[vertex]; n < new_num; n++) {
                        const auto move_idx = all_move_generator_indices_involving[vertex][n];
                        const auto& move = moves.get(move_idx);
                        this->set_active_in(move, move_idx, vertex);
                        vertices_in_updated_moves.insert(move.get_first_vertex());
                        vertices_in_updated_moves.insert(move.get_second_vertex());
                    }
                }

                curr_neighbors_num[vertex] = new_num;
            }

            // Note that before actually inserting the move generator indices in the list of move generators involving a vertex
            // we have to update the active/not active data structures that will be used in the following


            for (const auto vertex : vertices_to_update) {

                move_generator_indices_involving[vertex].clear();

                auto n = 0;
                for (; n < curr_neighbors_num[vertex]; n++) {
                    const auto idx = all_move_generator_indices_involving[vertex][n];
                    const auto& move = moves.get(idx);

                    const auto i = move.get_first_vertex();
                    const auto j = move.get_second_vertex();

                    assert(i == vertex || j == vertex);

                    assert(moves.get(idx + 1).get_first_vertex() == j);
                    assert(moves.get(idx + 1).get_second_vertex() == i);

                    assert(moves.get(idx).get_first_vertex() == i);
                    assert(moves.get(idx).get_second_vertex() == j);

                    assert(i < j);

                    if (vertex == i) {

                        move_generator_indices_involving[i].emplace_back(idx);

                        if (!is_active_in_other(move, idx, vertex)) {  // if it were not already active in j => j does not have it
                            move_generator_indices_involving[j].emplace_back(idx);
                        }

                    } else {

                        move_generator_indices_involving[j].emplace_back(idx);

                        if (!is_active_in_other(move, idx, vertex)) {  // if it were not already active in i => i does not have it
                            move_generator_indices_involving[i].emplace_back(idx);
                        }
                    }
                }

                for (; n < static_cast<int>(all_move_generator_indices_involving[vertex].size()); n++) {
                    const auto idx = all_move_generator_indices_involving[vertex][n];
                    const auto& move = moves.get(idx);

#ifndef NDEBUG
                    const auto i = move.get_first_vertex();
                    const auto j = move.get_second_vertex();
#endif

                    assert(i == vertex || j == vertex);

                    assert(moves.get(idx + 1).get_first_vertex() == j);
                    assert(moves.get(idx + 1).get_second_vertex() == i);

                    assert(moves.get(idx).get_first_vertex() == i);
                    assert(moves.get(idx).get_second_vertex() == j);

                    assert(i < j);

                    // add all move generators that are active in the other vertex
                    if (is_active_in_other(move, idx, vertex)) {
                        move_generator_indices_involving[vertex].emplace_back(idx);
                    }
                }
            }
        }

    public:
        KNeighborsMoveGeneratorsView(const cobra::Instance& inst, int k)
            : AbstractMoveGeneratorsView(inst,
                                         [&inst, k](auto i) -> std::vector<int> {
                                             auto endpoints = std::vector<int>();

                                             const auto max_neighbors = std::min(k, inst.get_vertices_num() - 1);

                                             for (auto n = 1, added = 0; added < max_neighbors; n++) {

                                                 assert(n < static_cast<int>(inst.get_neighbors_of(i).size()));

                                                 auto j = inst.get_neighbors_of(i)[n];

                                                 endpoints.emplace_back(j);
                                                 added++;
                                             }

                                             return endpoints;
                                         }),
              max_neighbors_num(std::min(k, inst.get_vertices_num() - 1))  // skip self-move
        {

            curr_neighbors_num.resize(inst.get_vertices_num(), 0);
        }
    };


}  // namespace cobra

#endif