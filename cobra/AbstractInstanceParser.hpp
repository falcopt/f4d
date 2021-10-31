#ifndef _F4D_ABSTRACTINSTANCEPARSER_HPP_
#define _F4D_ABSTRACTINSTANCEPARSER_HPP_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Flat2DVector.hpp"

namespace cobra {

    template <bool round_costs = true>
    class AbstractInstanceParser {

        Flat2DVector<float> costs_matrix;
        std::vector<std::vector<int>> neighbors;

    protected:
        const std::string& path;
        int vehicle_capacity;
        std::vector<float> x_coordinates;
        std::vector<float> y_coordinates;
        std::vector<int> demands;

    public:
        explicit AbstractInstanceParser(const std::string& path_) : path(path_) { }

        bool parse() {

            auto ok = parse_impl();

            if (!ok) {
                return ok;
            }

            auto matrix_size = static_cast<int>(demands.size());

            costs_matrix.resize(matrix_size, matrix_size);

            for (auto i = 0; i < matrix_size - 1; i++) {

                costs_matrix.at(i, i, 0.0f);

                for (auto j = i + 1; j < matrix_size; j++) {

                    costs_matrix.at(i, j,
                                    std::sqrt((x_coordinates[i] - x_coordinates[j]) * (x_coordinates[i] - x_coordinates[j]) +
                                              (y_coordinates[i] - y_coordinates[j]) * (y_coordinates[i] - y_coordinates[j])));

                    if constexpr (round_costs) {
                        costs_matrix.at(i, j, std::round(costs_matrix.at(i, j)));
                    }

                    costs_matrix.at(j, i, costs_matrix.at(i, j));
                }
            }

            neighbors.resize(matrix_size);

            auto all_vertices = std::vector<int>(matrix_size);

            for (auto i = 0; i < matrix_size; i++) {
                all_vertices[i] = i;
            }

            for (auto i = 0; i < matrix_size; i++) {

                std::sort(all_vertices.begin(), all_vertices.end(), [this, i](auto j, auto k) {
                    return costs_matrix.at(i, j) < costs_matrix.at(i, k);
                });

                neighbors[i] = all_vertices;

                // make sure the first vertex is i
                if (neighbors[i][0] != i) {
                    auto n = 1;
                    while (n < matrix_size) {
                        if (neighbors[i][n] == i) {
                            break;
                        }
                        n++;
                    }
                    std::swap(neighbors[i][0], neighbors[i][n]);
                }

                assert(neighbors[i][0] == i);
            }

            return ok;
        }

        virtual bool parse_impl() = 0;

        int get_vehicle_capacity() {
            return vehicle_capacity;
        }

        std::vector<float> get_x_coordinates() {
            return x_coordinates;
        }

        std::vector<float> get_y_coordinates() {
            return y_coordinates;
        }

        std::vector<int> get_demands() {
            return demands;
        }

        Flat2DVector<float> get_costs_matrix() {
            return costs_matrix;
        }

        std::vector<std::vector<int>> get_neighbors() {
            return neighbors;
        }
    };


    inline std::vector<std::string> split_line(const std::string& line, char separator) {
        auto tokens = std::vector<std::string>();

        auto stream = std::istringstream(line);
        auto token = std::string();

        while (std::getline(stream, token, separator)) {
            tokens.emplace_back(token);
        }

        return tokens;
    }

    template <bool round_costs>
    class TSPLIB95 : public AbstractInstanceParser<round_costs> {

    public:
        explicit TSPLIB95(const std::string& path_) : AbstractInstanceParser<round_costs>(path_) { }

        bool parse_impl() {

            try {

                auto stream = std::ifstream(this->path);

                if (!stream) {
                    return false;
                }

                auto line = std::string();

                auto matrix_size = 0;

                while (std::getline(stream, line)) {

                    if (line.rfind("DIMENSION") == 0) {

                        const auto tokens = split_line(line, ':');
                        matrix_size = std::stoi(tokens[1]);

                    } else if (line.rfind("CAPACITY") == 0) {

                        const auto tokens = split_line(line, ':');
                        this->vehicle_capacity = std::stoi(tokens[1]);

                    } else if (line.rfind("NODE_COORD_SECTION") == 0) {
                        break;
                    }
                }

                this->x_coordinates.resize(matrix_size);
                this->y_coordinates.resize(matrix_size);
                this->demands.resize(matrix_size);

                auto parse_depot_section = false;

                // coord section
                for (auto n = 0; n < matrix_size; n++) {

                    std::getline(stream, line);

                    if (line.rfind("DEMAND_SECTION") == 0) {
                        // depot was not here
                        parse_depot_section = true;

                        // move coordinates one position to the right

                        assert(n == matrix_size - 1);

                        for (auto i = matrix_size - 2; i >= 0; --i) {

                            this->x_coordinates[i + 1] = this->x_coordinates[i];
                            this->y_coordinates[i + 1] = this->y_coordinates[i];
                        }

                        goto demand_section;
                    }

                    std::istringstream iss(line);
                    float index, x, y;

                    iss >> index >> x >> y;

                    this->x_coordinates[n] = x;
                    this->y_coordinates[n] = y;
                }

                while (std::getline(stream, line)) {
                    if (line.rfind("DEMAND_SECTION") == 0) {
                        break;
                    }
                }

            demand_section:

                for (auto n = parse_depot_section ? 1 : 0; n < matrix_size; n++) {

                    std::getline(stream, line);

                    std::istringstream iss(line);
                    int index, q;

                    iss >> index >> q;

                    this->demands[n] = q;
                }

                if (parse_depot_section) {

                    while (std::getline(stream, line)) {
                        if (line.rfind("DEPOT_SECTION") == 0) {
                            break;
                        }
                    }

                    std::getline(stream, line);
                    std::istringstream iss(line);

                    float x, y;

                    iss >> x >> y;

                    this->x_coordinates[0] = x;
                    this->y_coordinates[0] = y;
                    this->demands[0] = 0;
                }

                return true;

            } catch (std::exception& e) {
                return false;
            }
            return true;
        }
    };


}  // namespace cobra


#endif