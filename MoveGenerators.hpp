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
        inline unsigned long get() { return value; }
        inline void increment() { ++value; }
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

        inline auto get_first_vertex() const { return first_vertex; }
        inline auto get_second_vertex() const { return second_vertex; }

        inline auto get_delta() const { return delta; }
        inline auto set_delta(float value) { delta = value; }

        inline auto get_heap_index() const { return heap_index; }
        inline auto set_heap_index(int index) { heap_index = index; }
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


        void reset() { BHeap::reset(); }
        bool is_empty() const { return BHeap::empty(); }
        void insert(MoveGenerator* mg) { BHeap::insert(mg); }
        MoveGenerator* get() { return BHeap::get(); }
        void remove(int heap_index) { BHeap::remove(heap_index); }
        void change_value(int heap_index, float value) { BHeap::update(heap_index, value); };
        int size() const { return BHeap::size(); };
        MoveGenerator* spy(int heap_index) { return BHeap::spy(heap_index); }

        static const int unheaped = -1;
    };


    /**
     * Move generators container
     */
    class AbstractMoveGeneratorsView;  // forward declaration
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
        MoveGenerators(const cobra::Instance& instance, std::vector<AbstractMoveGeneratorsView*>& views);

        virtual ~MoveGenerators();
        /**
         * Returns the move generator indexed `idx`
         * @param index
         * @return move generator indexed `idx`
         */
        inline MoveGenerator& get(int idx) { return moves[idx]; }
        /**
         * Set the percentage of move generators to consider in local search procedures
         * @param percentage float value in [0, 1]
         */
        void set_active_percentage(std::vector<float>& percentage, std::vector<int>& vertices);
        void set_active_percentage(float percentage);
        /**
         * Returns the heap used for storing move generators during local search applications
         * @return heap
         */
        inline MoveGeneratorsHeap& get_heap() { return heap; }
        /**
         * Returns the vector of move generators
         * @return vector of move generators
         */
        inline std::vector<MoveGenerator>& get_raw_vector() { return moves; }

        static inline auto get_twin_move_generator_index(int index) { return index ^ 1; }

        static inline auto get_base_move_generator_index(int index) { return index & ~1; }

        inline const auto& get_move_generator_indices_involving_1st(int vertex) const { return move_generator_indices_involving[vertex]; }

        inline auto get_move_generator_indices_involving_2nd(int vertex) const {
            const auto& v = move_generator_indices_involving[vertex];
            return VectorView<decltype(v.begin()), twin_functor>(v.begin(), v.end());
        }

        inline auto get_move_generator_indices_involving(int vertex) const {
            const auto& v = move_generator_indices_involving[vertex];
            return VectorView<decltype(v.begin()), base_functor>(v.begin(), v.end());
        }

        inline TimestampGenerator& get_timestamp_generator() { return timegen; }

        inline Flat2DVector<bool>& get_update_bits() { return update_bits; }

        inline std::vector<unsigned long>& get_vertex_timestamp() { return vertex_timestamp; }
    };

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

        void setup_active_trackers();

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
        AbstractMoveGeneratorsView(const cobra::Instance& instance_, std::function<std::vector<int>(int)> generator_);

        /**
         * Returns the generator
         * @return generator
         */
        auto get_generator() -> std::function<std::vector<int>(int)>&;
        auto get_move_generator_indices_involving(int vertex) -> std::vector<int>&;

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
                for (auto id : all_move_generator_indices_involving[i]) { set.insert(id); }
            }
            return 2 * static_cast<int>(set.size());  // here we are storing a single copy between (i, j) and (j,i)
        }

        inline auto is_active(int move_idx) {
            assert(move_map.count(move_idx));
            const auto mapped_idx = move_map[move_idx];
            return active_in[mapped_idx].first || active_in[mapped_idx].second;
        }
    };


    class KNeighborsMoveGeneratorsView : public AbstractMoveGeneratorsView {
        const int max_neighbors_num;          // k max, `set_active_percentage` can be used to vary 0 <= k <= k max
        std::vector<int> curr_neighbors_num;  // current k


        void __notify_move_indices(int vertex, std::vector<int>& indices, MoveGenerators& moves) override;
        void __notify_build_complete(MoveGenerators& a) override;
        void set_active_percentage(std::vector<float>& percentage, std::vector<int>& vertices, MoveGenerators& moves,
                                   cobra::VertexSet& vertices_in_updated_moves) override;

    public:
        KNeighborsMoveGeneratorsView(const cobra::Instance& inst, int k);
    };


    class CostBasedMoveGeneratorsView : public AbstractMoveGeneratorsView {

        std::vector<int> flat_indices;
        std::vector<float> curr_percentage;
        std::vector<std::vector<float>> inclusion_percentage;
        std::vector<int> curr_num;

        void __notify_move_indices([[maybe_unused]] int vertex, std::vector<int>& indices, MoveGenerators& moves) override;

        void __notify_build_complete(MoveGenerators& moves) override;

        void set_active_percentage(std::vector<float>& percentage, std::vector<int>& vertices, MoveGenerators& moves,
                                   cobra::VertexSet& vertices_in_updated_moves) override;

    public:
        CostBasedMoveGeneratorsView(const cobra::Instance& instance_, std::function<std::vector<int>(int)> generator_);
    };

    class _MoveGeneratorsHeap : private NonCopyable<_MoveGeneratorsHeap> {

        std::vector<MoveGenerator>& moves;
        int* heap;
        int heap_len;

        void heapify(int heap_index);
        void upsift(int heap_index);
        bool is_heap();
        void dump();

    public:
        explicit _MoveGeneratorsHeap(MoveGenerators& moves);
        virtual ~_MoveGeneratorsHeap();

        void reset();
        bool is_empty() const;
        void insert(int move_index);
        int get();
        void remove(int move_index);
        void change_value(int move_index, float value);
        int size() const;
        int spy(int heap_index) const;

        static const int unheaped;
    };

}  // namespace cobra

#endif