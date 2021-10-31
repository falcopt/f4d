#ifndef _F4D_NEIGHBORACCEPTANCE_HPP_
#define _F4D_NEIGHBORACCEPTANCE_HPP_

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

        template <bool maybe_infeasible>
        bool accept(const cobra::Solution &solution, const cobra::Solution &neighbor) {
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

        template <bool maybe_infeasible>
        bool accept(const cobra::Solution &solution, const cobra::Solution &neighbor) {
            return neighbor.get_cost() < solution.get_cost() + temperature;
        }

        float get_temperature() const {
            return temperature;
        }
    };


    class TimeBasedSimulatedAnnealing {

        float initial_temperature;
        float final_temperature;

        int max_time;

        std::mt19937 &rand_engine;
        std::uniform_real_distribution<float> uniform_dist;

        float a, b;


    public:
        TimeBasedSimulatedAnnealing(float initial_temperature_, float final_temperature_, std::mt19937 &rand_engine_, int max_time_)
            : rand_engine(rand_engine_), uniform_dist(0.0f, 1.0f) {

            initial_temperature = initial_temperature_;
            final_temperature = final_temperature_;
            max_time = max_time_;

            a = (initial_temperature - final_temperature) / (std::log(1.0f / static_cast<float>(max_time)));
            b = std::exp((-initial_temperature * std::log(static_cast<float>(max_time))) / (initial_temperature - final_temperature));

            // a = (final_temperature - initial_temperature) / static_cast<float>(max_time);
            // b = initial_temperature;
        }

        template <bool maybe_infeasible>
        bool accept(const cobra::Solution &solution, const cobra::Solution &neighbor, int elapsed_time) {
            const auto temperature = a * std::log(b * static_cast<float>(elapsed_time));
            // const auto temperature = a * static_cast<float>(elapsed_time) + b;
            return neighbor.get_cost() < solution.get_cost() - temperature * std::log(uniform_dist(rand_engine));
        }

        float get_temperature(int elapsed_time) const {
            return a * std::log(b * static_cast<float>(elapsed_time));
            // return a * elapsed_time + b;
        }
    };


    class TimeBasedDeterministicAnnealing {

        float initial_temperature;
        float final_temperature;
        int max_time;

        float a, b;

    public:
        TimeBasedDeterministicAnnealing(float initial_temperature_, float final_temperature_, int max_time_) {

            initial_temperature = initial_temperature_;
            final_temperature = final_temperature_;
            max_time = max_time_;

            a = (initial_temperature - final_temperature) / (std::log(1.0f / static_cast<float>(max_time)));
            b = std::exp((-initial_temperature * std::log(static_cast<float>(max_time))) / (initial_temperature - final_temperature));
        }

        template <bool maybe_infeasible>
        bool accept(const cobra::Solution &solution, const cobra::Solution &neighbor, int elapsed_time) {
            const auto temperature = a * std::log(b * static_cast<float>(elapsed_time));
            return neighbor.get_cost() < solution.get_cost() + temperature;
        }

        float get_temperature(int elapsed_time) const {
            return a * std::log(b * static_cast<float>(elapsed_time));
        }
    };


    class TimeBasedSimulatedAnnealing2 {

        float initial_temperature;
        float final_temperature;
        float temp_ratio;
        int period;

        std::mt19937 &rand_engine;
        std::uniform_real_distribution<float> uniform_dist;

        float factor;

    public:
        TimeBasedSimulatedAnnealing2(float initial_temperature_, float final_temperature_, std::mt19937 &rand_engine_, int max_iter)
            : rand_engine(rand_engine_), uniform_dist(0.0f, 1.0f) {

            initial_temperature = initial_temperature_;
            final_temperature = final_temperature_;
            temp_ratio = final_temperature / initial_temperature;
            period = max_iter;
        }


        template <bool maybe_infeasible>
        bool accept(const cobra::Solution &solution, const cobra::Solution &neighbor, double elapsed_time) {
            const double temperature = initial_temperature * std::pow(temp_ratio, elapsed_time / static_cast<float>(period));
            return neighbor.get_cost() < solution.get_cost() - temperature * std::log(uniform_dist(rand_engine));
        }

        float get_temperature(double elapsed_time) const {
            return initial_temperature * std::pow(temp_ratio, elapsed_time / static_cast<float>(period));
        }
    };

    class TimeBasedSimulatedAnnealing2Levels {

        double initial_temperature;
        double final_temperature;
        double initial_temperature_sect = 0.0;
        double temp_ratio;
        double temp_ratio_sect = 1.0;
        int period;
        int sector = -1;
        double sector_duration;

        std::mt19937 &rand_engine;
        std::uniform_real_distribution<double> uniform_dist;

        double factor;

    public:
        TimeBasedSimulatedAnnealing2Levels(double initial_temperature_, double final_temperature_, std::mt19937 &rand_engine_, int max_time,
                                           double sector_duration_)
            : rand_engine(rand_engine_), uniform_dist(0.0f, 1.0f) {

            initial_temperature = initial_temperature_;
            final_temperature = final_temperature_;
            temp_ratio = final_temperature / initial_temperature;
            period = max_time;
            sector_duration = sector_duration_;
        }


        template <bool maybe_infeasible>
        bool accept(const cobra::Solution &solution, const cobra::Solution &neighbor, double elapsed_time) {
            const auto temperature = get_temperature(elapsed_time);
            return neighbor.get_cost() < solution.get_cost() - temperature * std::log(uniform_dist(rand_engine));
        }

        double get_temperature(double elapsed_time) {
            int current_sector = static_cast<int>(elapsed_time / sector_duration);
            if (sector < current_sector) {
                sector = current_sector;
                initial_temperature_sect = initial_temperature * std::pow(temp_ratio, elapsed_time / static_cast<double>(period));
                temp_ratio_sect = final_temperature / initial_temperature_sect;
            }
            double elapsed_time_sect = elapsed_time - current_sector * sector_duration;
            const double temperature = initial_temperature_sect * std::pow(temp_ratio_sect, elapsed_time_sect / static_cast<double>(sector_duration));
            return temperature;
        }
    };


    class TimeBasedSimulatedAnnealing3Levels {

        double initial_temperature_1 = 0.0, initial_temperature_2 = 0.0, initial_temperature_3 = 0.0;
        double final_temperature;
        double temp_ratio_1 = 1.0, temp_ratio_2 = 1.0, temp_ratio_3 = 1.0;
        int sa1_tlimit, sa2_tlimit, sa3_tlimit;
        int sector_2 = -1, sector_3 = -1;
        double sector_duration;

        std::mt19937 &rand_engine;
        std::uniform_real_distribution<double> uniform_dist;

        double factor;

    public:
        TimeBasedSimulatedAnnealing3Levels(double initial_temperature_, double final_temperature_, std::mt19937 &rand_engine_, int max_time)
            : rand_engine(rand_engine_), uniform_dist(0.0f, 1.0f) {

            initial_temperature_1 = initial_temperature_;
            final_temperature = final_temperature_;
            temp_ratio_1 = 2 * final_temperature / initial_temperature_1;

            sa1_tlimit = max_time;
            sa2_tlimit = sa1_tlimit / 10;
            sa3_tlimit = sa2_tlimit / 10;
        }


        template <bool maybe_infeasible>
        bool accept(const cobra::Solution &solution, const cobra::Solution &neighbor, double elapsed_time) {
            const auto temperature = get_temperature(elapsed_time);
            return neighbor.get_cost() < solution.get_cost() - temperature * std::log(uniform_dist(rand_engine));
        }

        double get_temperature(double elapsed_time) {

            int current_sector_3 = elapsed_time / sa3_tlimit;
            if (sector_3 < current_sector_3) {
                sector_3 = current_sector_3;

                int current_sector_2 = elapsed_time / sa2_tlimit;
                if (sector_2 < current_sector_2) {
                    sector_2 = current_sector_2;
                    initial_temperature_2 = initial_temperature_1 * std::pow(temp_ratio_1, elapsed_time / static_cast<double>(sa1_tlimit));
                    temp_ratio_2 = final_temperature / initial_temperature_2;
                }
                double elapsed_time_sect_2 = elapsed_time - current_sector_2 * sa2_tlimit;
                initial_temperature_3 = initial_temperature_2 * std::pow(temp_ratio_2, elapsed_time_sect_2 / static_cast<double>(sa2_tlimit));
                temp_ratio_3 = final_temperature / initial_temperature_3;
            }
            double elapsed_time_sect_3 = elapsed_time - current_sector_3 * sa3_tlimit;
            double temperature = initial_temperature_3 * std::pow(temp_ratio_3, elapsed_time_sect_3 / static_cast<double>(sa3_tlimit));

            return temperature;
        }
    };


    class TimeBasedSimulatedAnnealingNLevels {

        double initial_temperature = 1.0;
        double final_temperature = 1.0;
        double temp_ratio = 1.0;
        std::vector<int> time_limits;
        std::vector<int> sect_counter;
        std::vector<double> init_temps;


        std::mt19937 &rand_engine;
        std::uniform_real_distribution<double> uniform_dist;

        double factor;

    public:
        TimeBasedSimulatedAnnealingNLevels(double initial_temperature_, double final_temperature_, std::mt19937 &rand_engine_, int max_time, int min_interval)
            : rand_engine(rand_engine_), uniform_dist(0.0f, 1.0f) {

            initial_temperature = initial_temperature_;
            final_temperature = final_temperature_;
            temp_ratio = final_temperature / initial_temperature;

            while (max_time > min_interval) {
                time_limits.push_back(max_time);
                max_time = (max_time + 1) / 2;
            }
            sect_counter = std::vector<int>(time_limits.size(), -1);
            init_temps = std::vector<double>(time_limits.size(), initial_temperature);
        }


        template <bool maybe_infeasible>
        bool accept(const cobra::Solution &solution, const cobra::Solution &neighbor, double elapsed_time) {
            const auto temperature = get_temperature(elapsed_time);
            return neighbor.get_cost() < solution.get_cost() - temperature * std::log(uniform_dist(rand_engine));
        }

        double get_temperature(double elapsed_time) {
            return get_temperature(elapsed_time, time_limits.size() - 1);
        }

        double get_temperature(double elapsed_time, size_t sect_i) {
            int current_sector = elapsed_time / time_limits[sect_i];
            if (sect_i > 0 && sect_counter[sect_i] < current_sector) {
                sect_counter[sect_i] = current_sector;
                init_temps[sect_i] = get_temperature(elapsed_time, sect_i - 1);
            }
            double elapsed_time_sect = elapsed_time - current_sector * time_limits[sect_i];
            const auto factor = final_temperature / init_temps[sect_i];
            const auto temperature = init_temps[sect_i] * std::pow(factor, elapsed_time_sect / static_cast<double>(time_limits[sect_i]));
            if (sect_i < time_limits.size() - 1) {
                std::cout << "temperature" << temperature << "= init_temps[sect_i]" << init_temps[sect_i] << " * std::pow(factor" << factor
                          << ", elapsed_time_sect" << elapsed_time_sect << " / static_cast<double>(time_limits[sect_i]" << time_limits[sect_i] << "));\n";
                std::cout << sect_i << ": InitialTemperature:" << init_temps[sect_i] << ", Temperature: " << temperature << '\n';
            }
            return temperature;
        }
    };


    class TimeBasedSimulatedAnnealinglastImprovement {

        double initial_temperature;
        double final_temperature;
        double initial_temperature_sect = 0.0;
        double temp_ratio;
        double temp_ratio_sect = 1.0;
        int period;
        int sector = -1;
        double sector_duration;
        double last_temp;
        double temperature;

        std::mt19937 &rand_engine;
        std::uniform_real_distribution<double> uniform_dist;

        double factor;

    public:
        TimeBasedSimulatedAnnealinglastImprovement(double initial_temperature_, double final_temperature_, std::mt19937 &rand_engine_, int max_time,
                                                   double sector_duration_)
            : rand_engine(rand_engine_), uniform_dist(0.0f, 1.0f) {

            initial_temperature = initial_temperature_;
            last_temp = initial_temperature;
            initial_temperature_sect = initial_temperature;
            final_temperature = final_temperature_;
            temp_ratio = final_temperature / initial_temperature;
            period = max_time;
            sector_duration = sector_duration_;
        }


        template <bool maybe_infeasible>
        bool accept(const cobra::Solution &solution, const cobra::Solution &neighbor, double elapsed_time) {
            const auto temperature = get_temperature(elapsed_time);
            return neighbor.get_cost() < solution.get_cost() - temperature * std::log(uniform_dist(rand_engine));
        }

        double get_temperature(double elapsed_time) {
            int current_sector = static_cast<int>(elapsed_time / sector_duration);
            if (sector < current_sector) {
                sector = current_sector;
                initial_temperature_sect = 20 * last_temp;
                std::cerr << "Temps: " << temperature << " LastImpr: " << last_temp << " Init: " << initial_temperature_sect << "\n";
                temp_ratio_sect = final_temperature / initial_temperature_sect;
            }
            double elapsed_time_sect = elapsed_time - current_sector * sector_duration;
            temperature = initial_temperature_sect * std::pow(temp_ratio_sect, elapsed_time_sect / static_cast<double>(sector_duration));
            return temperature;
        }

        void save_last_improving_temp() {
            last_temp = temperature;
        }
    };


}  // namespace cobra

#endif  // COBRA_INCLUDE_COBRA_NEIGHBORACCEPTANCE_HPP_