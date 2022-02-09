#ifndef _F4D_NEIGHBORACCEPTANCE_HPP_
#define _F4D_NEIGHBORACCEPTANCE_HPP_

#include <cmath>
#include <random>

#include "Solution.hpp"

namespace cobra {

    class SimulatedAnnealing {

        float initial_temperature;
        float final_temperature;
        float temperature;
        int period;

        std::mt19937 &rand_engine;
        std::uniform_real_distribution<float> uniform_dist;

        float factor;

    public:
        SimulatedAnnealing(float initial_temperature_, float final_temperature_, std::mt19937 &rand_engine_, int max_iter)
            : rand_engine(rand_engine_), uniform_dist(0.0f, 1.0f) {

            initial_temperature = initial_temperature_;
            final_temperature = final_temperature_;
            period = max_iter;

            temperature = initial_temperature;
            factor = std::pow(final_temperature / initial_temperature, 1.0f / static_cast<float>(period));
        }

        void decrease_temperature() {
            temperature *= factor;
        }

        bool accept(const Solution &solution, const Solution &neighbor) {
            return neighbor.get_cost() < solution.get_cost() - temperature * std::log(uniform_dist(rand_engine));
        }

        float get_temperature() const {
            return temperature;
        }
    };


    class DeterministicAnnealing {

        float initial_temperature;
        float final_temperature;
        float temperature;
        int period;

        float factor;

    public:
        DeterministicAnnealing(float initial_temperature_, float final_temperature_, int max_iter) {

            initial_temperature = initial_temperature_;
            final_temperature = final_temperature_;
            period = max_iter;

            temperature = initial_temperature;
            factor = std::pow(final_temperature / initial_temperature, 1.0f / static_cast<float>(period));
        }

        void decrease_temperature() {
            temperature *= factor;
        }

        bool accept(const Solution &solution, const Solution &neighbor) {
            return neighbor.get_cost() < solution.get_cost() + temperature;
        }

        float get_temperature() const {
            return temperature;
        }
    };

    class TimeBasedSimulatedAnnealing {

        float initial_temperature;
        float final_temperature;
        float temp_ratio;
        int period;

        std::mt19937 &rand_engine;
        std::uniform_real_distribution<float> uniform_dist;

        float factor;

    public:
        TimeBasedSimulatedAnnealing(float initial_temperature_, float final_temperature_, std::mt19937 &rand_engine_, int max_iter)
            : rand_engine(rand_engine_), uniform_dist(0.0f, 1.0f) {

            initial_temperature = initial_temperature_;
            final_temperature = final_temperature_;
            temp_ratio = final_temperature / initial_temperature;
            period = max_iter;
        }


        bool accept(const Solution &solution, const Solution &neighbor, double elapsed_time) {
            const double temperature = initial_temperature * std::pow(temp_ratio, elapsed_time / static_cast<float>(period));
            return neighbor.get_cost() < solution.get_cost() - temperature * std::log(uniform_dist(rand_engine));
        }

        float get_temperature(double elapsed_time) const {
            return initial_temperature * std::pow(temp_ratio, elapsed_time / static_cast<float>(period));
        }
    };


}  // namespace cobra

#endif