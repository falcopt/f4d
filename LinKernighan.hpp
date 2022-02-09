#ifndef _F4D_LINKERNIGHAN_HPP_
#define _F4D_LINKERNIGHAN_HPP_

#include <queue>
#include <set>

#include "cobra/Instance.hpp"
#include "cobra/Solution.hpp"

#define MIN_TOUR_SIZE (25)
#define MAX_NEIGHBORS_IN_TOUR (10)

template <bool maybe_infeasible = false>
class LinKernighan {

    struct Flip {
        int begin;
        int end;
        Flip(int begin_, int end_) : begin(begin_), end(end_) { }
    };


    class Tour : NonCopyable<Tour> {

    public:
        static const bool left_to_right = false;
        static const bool right_to_left = true;

        explicit Tour(const cobra::Instance& instance_) : instance(instance_) {

            tour.reserve(instance.get_vertices_num());
            position.resize(instance.get_vertices_num());
            neighbors.resize(instance.get_vertices_num());
            for (auto i = instance.get_vertices_begin(); i < instance.get_vertices_end(); i++) {
                neighbors[i].resize(MAX_NEIGHBORS_IN_TOUR);
            }

            cost = 0.0;
            orientation = Tour::left_to_right;
            tour_size = 0;
        }

        void build_from(const cobra::Solution& solution, const int route) {

            tour.clear();

            cost = solution.get_route_cost(route);
            orientation = Tour::left_to_right;
            tour_size = solution.get_route_size(route) + 1;

            tour.emplace_back(instance.get_depot());
            position[instance.get_depot()] = 0;

            int index = 1;

            for (auto i = solution.get_first_customer(route); i != instance.get_depot(); i = solution.get_next_vertex(i)) {
                tour.emplace_back(i);
                position[i] = index++;
            }

            auto flags = std::vector<bool>(tour.size());

            for (auto i : tour) {
                neighbors[i].clear();
                std::fill(flags.begin(), flags.end(), false);
                flags[position[i]] = true;
                for (auto m = 0; m < MAX_NEIGHBORS_IN_TOUR; m++) {
                    auto min_index = -1;
                    auto min_value = std::numeric_limits<float>::max();
                    for (auto n = 0u; n < tour.size(); n++) {
                        if (flags[n]) {
                            continue;
                        }
                        const auto vertex = tour[n];
                        if (instance.get_cost(i, vertex) < min_value) {
                            min_value = instance.get_cost(i, vertex);
                            min_index = n;
                        }
                    }
                    flags[min_index] = true;
                    neighbors[i].emplace_back(tour[min_index]);
                }
                // std::sort(neighbors[i].begin(), neighbors[i].end(), [this, i](auto a, auto b){return instance.get_cost(a, i) < instance.get_cost(b, i); });
            }
        }

        /**
         * Inverts the segment of the tour from a to b possibly changing the tour orientation
         */
        float flip(int a, int b) {

            const auto direct_path_length = vertices_in_path(a, b);
            const auto inverse_path_length = tour_size - direct_path_length;

            const auto prev_a = prev(a);
            const auto next_b = next(b);

            int left;
            int right;
            int swaps;

            // check whether it is more convenient to flip a-b or next_b-prev_a
            if (direct_path_length <= inverse_path_length) {

                left = a;
                right = b;
                swaps = direct_path_length / 2;

            } else {

                left = prev_a;
                right = next_b;
                swaps = inverse_path_length / 2;

                orientation = !orientation;  // toggle orientation
            }

            if (!swaps) {
                return 0;
            }

            while (swaps--) {

                const auto next_left = next(left);
                const auto prev_right = prev(right);

                const auto pleft = position[left];
                const auto pright = position[right];

                tour[pleft] = right;
                tour[pright] = left;
                position[left] = pright;
                position[right] = pleft;

                left = next_left;
                right = prev_right;
            }

            const auto delta = -instance.get_cost(prev_a, a) - instance.get_cost(b, next_b) + instance.get_cost(prev_a, b) + instance.get_cost(next_b, a);

            cost += delta;


            return delta;
        }

        /**
         * Returns the vertex immediately after a in the tour
         * @param a vertex
         * @return next(a)
         */
        int next(int a) const {
            if (!orientation) {  // orientation == Tour::left_to_right
                return tour[(position[a] + 1) % tour_size];
            } else {
                return tour[((position[a] - 1) % tour_size) + ((position[a] >= 1) ? 0 : tour_size)];
            }
        }

        /**
         * Returns the vertex immediately before a in the tour
         * @param a vertex
         * @return prev(a)
         */
        int prev(int a) const {
            if (orientation) {  // orientation == Tour::right_to_left
                return tour[(position[a] + 1) % tour_size];
            } else {
                return tour[((position[a] - 1) % tour_size) + ((position[a] >= 1) ? 0 : tour_size)];
            }
        }

        /**
         * Returns true if b lies in the segment of the tour from a to c, and return false otherwise
         * @param a left vertex
         * @param b vertex
         * @param c right vertex
         * @return b is between a and c
         */
        bool sequence(int a, int b, int c) const {

            if (orientation) {  // orientation == Tour::right_to_left
                std::swap(a, c);
            }

            const auto pa = position[a];
            const auto pb = position[b];
            const auto pc = position[c];
            if (pa <= pc) {
                return pb >= pa && pb <= pc;
            } else {
                return pb >= pa || pb <= pc;
            }
        }

        /**
         * Forces a specific tour orientation
         * @param orientation_
         */
        void set_orientation(bool orientation_) {
            orientation = orientation_;
        }

        /**
         * Computes the distance in terms of vertices between a and b
         * @param a vertex
         * @param b vertex
         * @return distance
         */
        int vertices_in_path(int a, int b) const {

            if (orientation) {  // orientation == Tour::right_to_left
                std::swap(a, b);
            }

            const auto pa = position[a];
            const auto pb = position[b];

            const auto diff = pb - pa;

            if (diff >= 0) {
                return 1 + diff;
            } else {
                return 1 + pb + (tour_size - pa);
            }
        }

        /**
         * Checks tour feasibility (internal debug procedure)
         * @return whether the tour is feasible
         */
        bool is_feasible() {

            auto flags = std::vector<bool>(instance.get_vertices_num());

            // Ensure vertices are visited once and next and prev work correctly
            std::fill(flags.begin(), flags.end(), false);
            for (auto curr = instance.get_depot(), n = 0; n < tour_size; n++, curr = next(curr)) {
                if (flags[curr]) {
                    std::cout << "Vertex " << curr << " visited more than once while traversing with NEXT\n";
                    return false;
                }
                flags[curr] = true;
            }

            std::fill(flags.begin(), flags.end(), false);
            for (auto curr = instance.get_depot(), n = 0; n < tour_size; n++, curr = prev(curr)) {
                if (flags[curr]) {
                    std::cout << "Vertex " << curr << " visited more than once while traversing with PREV\n";
                    return false;
                }
                flags[curr] = true;
            }

            auto tour_load = 0;

            for (auto curr = instance.get_depot(), n = 0; n < tour_size; n++, curr = next(curr)) {
                tour_load += instance.get_demand(curr);
            }

            // check the stored tour cost
            auto recomputed = 0.0;
            for (auto curr = instance.get_depot(), n = 0; n < tour_size; n++, curr = next(curr)) {
                recomputed += instance.get_cost(curr, next(curr));
            }

            if (std::fabs(cost - recomputed) > 0.01) {
                std::cout << "Stored cost " << cost << " is different from computed one " << recomputed << "\n";
                return false;
            }

            return true;
        }

        float get_cost() {
            return static_cast<float>(cost);
        }

        void print() {
            for (auto curr = instance.get_depot(), n = 0; n < tour_size; n++, curr = next(curr)) {
                std::cout << curr << " ";
            }
            std::cout << "\n";
        }

        const std::vector<int>& get_tour() {
            return tour;
        }

        const std::vector<int>& get_neighbors(int vertex) {
            return neighbors[vertex];
        }

    private:
        const cobra::Instance& instance;
        double cost;
        int tour_size;  // DO NOT USE UNSIGNED INT!

        std::vector<int> tour;
        std::vector<int> position;
        std::vector<std::vector<int>> neighbors;

        bool orientation;
    };

    const cobra::Instance& instance;

    const float tolerance;

    Tour tour;

    void lin_kernighan() {

        // mark all vertices
        auto marked = std::vector<bool>(instance.get_vertices_num(), true);
        auto bases = std::queue<int>();
        for (auto i : tour.get_tour()) {
            bases.push(i);
        }

        // while there exist marked vertices
        while (!bases.empty()) {

            // select a marked vertex
            const auto base = bases.front();

            assert(tour.is_feasible());

            // auto old_cost = tour.get_cost();

            // call lk_search
            const auto flips = lk_search(base);

            assert(tour.is_feasible());

            // if an improving sequence of flips is found
            if (!flips.empty()) {

                for (const auto flip : flips) {

                    if (!marked[flip.begin]) {
                        marked[flip.begin] = true;
                        bases.push(flip.begin);
                    }

                    if (!marked[flip.end]) {
                        marked[flip.end] = true;
                        bases.push(flip.end);
                    }
                }


            } else {

                // unmark v
                bases.pop();
            }

            // assert(tour.get_cost() <= old_cost);
        }
    }

    std::vector<Flip> lk_search(int base) {

        auto flips = std::vector<Flip>();

        step(base, 1, 0.0f, tour.get_cost(), flips);

        if (!flips.empty()) {
            return flips;
        }

        alternate_step(base, tour.get_cost(), flips);

        assert(tour.is_feasible());

        return flips;
    }

    void step(int base, int level, float delta, float initial_cost, std::vector<Flip>& flips) {

        enum MoveType { LK, MAK_MORTON };

        struct CandidateMove {
            int endpoint;
            float value;
            MoveType type;
            CandidateMove(int endpoint_, float value_, MoveType type_) : endpoint(endpoint_), value(value_), type(type_) { }
        };

        // create the lk-ordering for base
        auto lk_ordering = std::vector<CandidateMove>();
        for (auto a : tour.get_neighbors(tour.next(base))) {
            const auto probe = tour.prev(a);
            if (probe == base || probe == tour.next(base) || probe == tour.prev(base)) {
                continue;
            }
            if (delta + instance.get_cost(base, tour.next(base)) - instance.get_cost(tour.next(base), a) <= 0) {
                continue;
            }
            // probe is a promising vertex for a LK move
            lk_ordering.emplace_back(a, instance.get_cost(probe, a) - instance.get_cost(tour.next(base), a), MoveType::LK);
        }

        for (auto a : tour.get_neighbors(base)) {
            if (a == tour.next(base) || a == tour.prev(base)) {
                continue;
            }
            if (delta + instance.get_cost(base, tour.next(base)) - instance.get_cost(base, a) <= 0) {
                continue;
            }
            lk_ordering.emplace_back(a, instance.get_cost(a, tour.next(a)) - instance.get_cost(base, a), MoveType::MAK_MORTON);
        }

        std::sort(lk_ordering.begin(), lk_ordering.end(), [](const auto& a, const auto& b) {
            return a.value > b.value;
        });

        auto i = 1;

        const auto breadth = [](int lev) {
            if (lev <= 2) {
                return 5;
            } else {
                return 1;
            }
        };

        for (auto n = 0u; n < lk_ordering.size() && i < breadth(level); n++) {

            const auto a = lk_ordering[n].endpoint;

            if (lk_ordering[n].type == MoveType::LK) {

                flips.emplace_back(tour.next(base), tour.prev(a));

                const auto g = tour.flip(tour.next(base), tour.prev(a));

                step(base, level + 1, delta + g, initial_cost, flips);

            } else {

                const auto new_base = tour.next(a);
                const auto old_base = base;

                flips.emplace_back(new_base, base);

                const auto g = tour.flip(new_base, base);

                base = new_base;

                step(base, level + 1, delta + g, initial_cost, flips);

                base = old_base;
            }


            if (tour.get_cost() + tolerance < initial_cost) {
                return;
            } else {

                const auto flip = flips.back();
                tour.flip(flip.end, flip.begin);
                flips.pop_back();
                i++;
            }
        }
    }

    void alternate_step(int base, float initial_cost, std::vector<Flip>& flips) {

        const auto s1 = tour.next(base);

        struct CandidateVertex {
            int vertex;
            float value;
            CandidateVertex(int vertex_, float value_) : vertex(vertex_), value(value_) { }
        };

        struct CandidateEdge {
            int vertex1;
            int vertex2;
            float value;
            CandidateEdge(int vertex1_, int vertex2_, float value_) : vertex1(vertex1_), vertex2(vertex2_), value(value_) { }
        };

        // create the A-ordering from the neighbors of s1
        auto a_ordering = std::vector<CandidateVertex>();
        for (auto a : tour.get_neighbors(s1)) {
            const auto a1 = tour.next(a);
            // in selecting a, we consider only the promising neighbors of s1
            // if(probe == base || probe == tour.next(base) || probe == tour.prev(base)) { continue; } TODO qualche check?
            // if(a1 == base || a1 == s1 || a1 == tour.prev(base)) { continue; }
            if (instance.get_cost(base, s1) - instance.get_cost(s1, a) <= 0) {
                continue;
            }
            a_ordering.emplace_back(a, instance.get_cost(a1, a) - instance.get_cost(s1, a));
        }
        std::sort(a_ordering.begin(), a_ordering.end(), [](const auto& a, const auto& b) {
            return a.value > b.value;
        });

        auto i = 1;

        static const auto a_breadth = 10;
        for (auto a_n = 0u; a_n < a_ordering.size() && i <= a_breadth; a_n++) {

            const auto a = a_ordering[a_n].vertex;
            const auto a1 = tour.next(a);

            // create the B-ordering from the neighbors of a1
            auto b_ordering = std::vector<CandidateEdge>();
            for (auto b : tour.get_neighbors(a1)) {
                if (b == base || b == s1 || b == a) {
                    continue;
                }
                if (instance.get_cost(tour.next(a), b) >= instance.get_cost(a, a1) + instance.get_cost(base, s1) - instance.get_cost(s1, a)) {
                    continue;
                }
                const auto b1 = tour.prev(b);
                b_ordering.emplace_back(b, b1, instance.get_cost(b1, b) - instance.get_cost(a, b));
            }
            std::sort(b_ordering.begin(), b_ordering.end(), [](const auto& x, const auto& y) {
                return x.value > y.value;
            });

            auto j = 1;

            static const auto b_breadth = 10;
            for (auto b_n = 0u; b_n < b_ordering.size() && j <= b_breadth; b_n++) {

                const auto b = b_ordering[b_n].vertex1;
                const auto b1 = b_ordering[b_n].vertex2;

                if (tour.sequence(s1, b, a)) {

                    flips.emplace_back(s1, b);
                    flips.emplace_back(b, a);
                    flips.emplace_back(s1, a);
                    flips.emplace_back(b, s1);
                    flips.emplace_back(a, b1);

                    const auto delta = tour.flip(s1, b) + tour.flip(b, a) + tour.flip(s1, a) + tour.flip(b, s1) + tour.flip(a, b1);

                    assert(tour.is_feasible());

                    step(base, 3, delta, initial_cost, flips);

                    if (tour.get_cost() + tolerance < initial_cost) {

                        return;

                    } else {

                        tour.flip(b1, a);
                        tour.flip(s1, b);
                        tour.flip(a, s1);
                        tour.flip(a, b);
                        tour.flip(b, s1);

                        assert(tour.is_feasible());

                        flips.pop_back();
                        flips.pop_back();
                        flips.pop_back();
                        flips.pop_back();
                        flips.pop_back();
                    }

                } else {

                    // create the D-ordering from the neighbors of b1
                    auto d_ordering = std::vector<CandidateEdge>();
                    for (auto d : tour.get_neighbors(b1)) {
                        // d must lie on the segment s1-a
                        if (!tour.sequence(s1, d, a)) {
                            continue;
                        }
                        if (d == base || d == s1 || d == a || d == a1 || d == b) {
                            continue;
                        }
                        if (instance.get_cost(b1, d) >= instance.get_cost(b, b1) + instance.get_cost(base, s1) - instance.get_cost(s1, a) +
                                                            instance.get_cost(a, a1) - instance.get_cost(a1, b)) {
                            continue;
                        }
                        const auto d1 = tour.prev(d);
                        d_ordering.emplace_back(d, d1, instance.get_cost(d1, d) - instance.get_cost(b1, d));
                    }
                    std::sort(d_ordering.begin(), d_ordering.end(), [](const auto& x, const auto& y) {
                        return x.value > y.value;
                    });

                    auto k = 1;

                    static const auto d_breadth = 10;
                    for (auto d_n = 0u; d_n < d_ordering.size() && k <= d_breadth; d_n++) {

                        const auto d = d_ordering[d_n].vertex1;
                        const auto d1 = d_ordering[d_n].vertex2;

                        flips.emplace_back(s1, b1);
                        flips.emplace_back(b1, d1);
                        flips.emplace_back(a1, s1);
                        flips.emplace_back(s1, d1);
                        flips.emplace_back(d, a);
                        flips.emplace_back(a1, b1);

                        const auto delta = tour.flip(s1, b1) + tour.flip(b1, d1) + tour.flip(a1, s1) + tour.flip(s1, d1) + tour.flip(d, a) + tour.flip(a1, b1);

                        assert(tour.is_feasible());

                        step(base, 4, delta, initial_cost, flips);

                        if (tour.get_cost() + tolerance < initial_cost) {

                            return;

                        } else {

                            tour.flip(b1, a1);
                            tour.flip(a, d);
                            tour.flip(d1, s1);
                            tour.flip(s1, a1);
                            tour.flip(d1, b1);
                            tour.flip(b1, s1);

                            assert(tour.is_feasible());

                            flips.pop_back();
                            flips.pop_back();
                            flips.pop_back();
                            flips.pop_back();
                            flips.pop_back();
                            flips.pop_back();

                            k++;
                        }
                    }
                }

                j++;
            }

            i++;
        }
    }

    void apply(cobra::Solution& solution, int route) {

        tour.build_from(solution, route);

        const auto initial_cost = tour.get_cost();

        lin_kernighan();

        if (tour.get_cost() + tolerance < initial_cost) {

            for (auto vertex : tour.get_tour()) {
                if (vertex == this->instance.get_depot()) {
                    continue;
                }
                solution.remove_vertex(route, vertex);
            }

            const auto first_customer = tour.next(this->instance.get_depot());
            for (auto curr = (first_customer); curr != this->instance.get_depot(); curr = tour.next(curr)) {
                solution.insert_vertex_before(route, this->instance.get_depot(), curr);
            }

            assert(solution.is_feasible());
            assert(std::fabs(solution.get_route_cost(route) - tour.get_cost()) < 0.01f);
        }
    }

public:
    LinKernighan(const cobra::Instance& instance_, const float tolerance_) : instance(instance_), tolerance(tolerance_), tour(instance) { }

    void apply(cobra::Solution& solution) {

        auto routes = std::set<int>();

        for (auto i = solution.get_cache_begin(); i != solution.get_cache_end(); i = solution.get_cache_next(i)) {
            if (i == instance.get_depot()) {
                continue;
            }
            routes.insert(solution.get_route_index(i));
        }

        for (auto route : routes) {
            if (solution.get_route_size(route) >= MIN_TOUR_SIZE) {  // do not apply LK on short routes ..
                apply(solution, route);
            }
        }
    }
};


#endif