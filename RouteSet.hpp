#ifndef _F4D_ROUTESET_HPP_
#define _F4D_ROUTESET_HPP_

#include "Instance.hpp"
#include "Solution.hpp"
namespace cobra {

#define SET_PARTITIONING_MAX_POOL_SIZE (100000)

    // Store some customers as a set with a cost label associated
    template <bool maybe_infeasible = false>
    class RouteSet {

        const cobra::Instance &instance;

        unsigned int words_num = 0;
        unsigned int *word = nullptr;
        unsigned int *bit = nullptr;

        static const unsigned int HASH_LOAD = 257;
        static const unsigned int BUCKETS_NUM = 7057;
        unsigned int max_route_len = 0;

        unsigned int (*head)[HASH_LOAD][BUCKETS_NUM] = nullptr;
        unsigned int *next = nullptr;

        unsigned short *route_load = nullptr;
        unsigned int **route_set = nullptr;

        unsigned int entries_num = 0;

    public:
        explicit RouteSet(const cobra::Instance &instance_) : instance(instance_) {

            words_num = static_cast<unsigned int>(std::ceil(instance.get_vertices_num() / 30.0));
            word = cobra::request_raw_contiguous_memory<unsigned int>(instance.get_vertices_num());
            bit = cobra::request_raw_contiguous_memory<unsigned int>(instance.get_vertices_num());

            auto pow = 0;
            for (auto i = instance.get_vertices_begin(); i < instance.get_vertices_end(); i++) {
                word[i] = static_cast<unsigned int>(i / 30.0);
                bit[i] = static_cast<unsigned int>(std::pow(2.0, pow));
                pow++;
                if (pow == 30) {
                    pow = 0;
                }
            }

            auto demands = std::vector<int>();
            demands.reserve(instance.get_customers_num());

            for (auto i = instance.get_customers_begin(); i < instance.get_customers_end(); i++) {
                demands.emplace_back(instance.get_demand(i));
            }

            std::sort(demands.begin(), demands.end(), [](const short i, const short j) {
                return i < j;
            });

            auto sum = 0;
            for (max_route_len = 0; max_route_len < static_cast<unsigned int>(instance.get_vertices_num()); max_route_len++) {
                sum += demands[max_route_len];
                if (sum > instance.get_vehicle_capacity()) {
                    break;
                }
            }

            head = cobra::request_raw_contiguous_memory<unsigned int[HASH_LOAD][BUCKETS_NUM]>(instance.get_vertices_num());
            next = cobra::request_raw_contiguous_memory<unsigned int>(SET_PARTITIONING_MAX_POOL_SIZE);

            entries_num = 0;

            route_seq = cobra::request_raw_contiguous_memory<unsigned short>(SET_PARTITIONING_MAX_POOL_SIZE, instance.get_vertices_num());
            route_load = cobra::request_raw_contiguous_memory<unsigned short>(SET_PARTITIONING_MAX_POOL_SIZE);
            route_cost = cobra::request_raw_contiguous_memory<float>(SET_PARTITIONING_MAX_POOL_SIZE);
            route_len = cobra::request_raw_contiguous_memory<unsigned short>(SET_PARTITIONING_MAX_POOL_SIZE);
            route_set = cobra::request_raw_contiguous_memory<unsigned int>(SET_PARTITIONING_MAX_POOL_SIZE, words_num);
            solution_cost = cobra::request_raw_contiguous_memory<float>(SET_PARTITIONING_MAX_POOL_SIZE);
        }

        void clear() {
            entries_num = 0;

            for (auto i = 0; i < SET_PARTITIONING_MAX_POOL_SIZE; i++) {
                for (auto j = 0u; j < words_num; j++) {
                    route_set[i][j] = 0;
                }
            }

            for (auto i = 0u; i < HASH_LOAD; i++) {
                for (auto j = 0u; j < BUCKETS_NUM; j++) {
                    for (auto k = 0; k < instance.get_vertices_num(); k++) {
                        head[i][j][k] = 0;
                    }
                }
            }
        }

        virtual ~RouteSet() {

            cobra::release_raw_contiguous_memory(word);
            cobra::release_raw_contiguous_memory(bit);

            cobra::release_raw_contiguous_memory(head);
            cobra::release_raw_contiguous_memory(next);

            cobra::release_raw_contiguous_memory(route_seq);
            cobra::release_raw_contiguous_memory(route_load);
            cobra::release_raw_contiguous_memory(route_cost);
            cobra::release_raw_contiguous_memory(route_len);
            cobra::release_raw_contiguous_memory(route_set);
            cobra::release_raw_contiguous_memory(solution_cost);
        }

        inline unsigned int size() const {
            return entries_num;
        }

        inline unsigned int begin() const {
            return 0;
        }

        inline unsigned int end() const {
            return entries_num;
        }

        float *route_cost = nullptr;
        unsigned short **route_seq = nullptr;
        unsigned short *route_len = nullptr;
        float *solution_cost = nullptr;  // cost of the best solution in which the route was found

        void add_routes(cobra::Solution<maybe_infeasible> &solution) {
            for (auto route = solution.get_first_route(); route != Solution<maybe_infeasible>::dummy_route; route = solution.get_next_route(route)) {
                lookup_and_insert_if_not_exists(solution, route);
            }
        }

        unsigned int lookup_and_insert_if_not_exists(cobra::Solution<maybe_infeasible> &solution, int route) {

            if (entries_num == SET_PARTITIONING_MAX_POOL_SIZE) {
                std::cout << "Reached max pool size. Skipping route.\n";
                return 0;
            }

            assert(entries_num + 1 < SET_PARTITIONING_MAX_POOL_SIZE);

            const auto insertion_position = entries_num;

            route_len[insertion_position] = 0;
            route_cost[insertion_position] = solution.get_route_cost(route);
            route_load[insertion_position] = solution.get_route_load(route);
            solution_cost[insertion_position] = solution.get_cost();

            for (auto customer = solution.get_first_customer(route); customer != instance.get_depot(); customer = solution.get_next_vertex(customer)) {
                route_set[insertion_position][word[customer]] |= bit[customer];
                route_seq[insertion_position][route_len[insertion_position]++] = static_cast<unsigned short>(customer);
            }

            auto key = 0ul;
            for (auto w = 0u; w < words_num; w++) {
                key += route_set[insertion_position][w];
            }
            key %= BUCKETS_NUM;

            const auto l = solution.get_route_load(route) % HASH_LOAD;
            const auto j = solution.get_route_size(route) % max_route_len;

            auto curr = head[j][l][key];

            while (curr) {

                if (route_load[curr] != route_load[insertion_position] || route_len[curr] != route_len[insertion_position]) {
                    goto notFound;
                }

                for (auto w = 0u; w < words_num; w++) {
                    if (route_set[curr][w] != route_set[insertion_position][w]) {
                        goto notFound;
                    }
                }

                //  Found!

                for (auto w = 0u; w < words_num; w++) {
                    route_set[insertion_position][w] = 0;
                }

                if (route_cost[insertion_position] + 0.01f < route_cost[curr]) {

                    route_cost[curr] = route_cost[insertion_position];
                    for (auto n = 0; n < route_len[curr]; n++) {
                        route_seq[curr][n] = route_seq[insertion_position][n];
                    }
                }

                if (solution_cost[insertion_position] + 0.01f < solution_cost[curr]) {
                    solution_cost[curr] = solution_cost[insertion_position];
                }

                return curr;

            notFound:
                curr = next[curr];
            }

            // not found
            next[insertion_position] = head[j][l][key];
            head[j][l][key] = insertion_position;

            entries_num++;

            return insertion_position;
        }
    };

}  // namespace cobra

#endif