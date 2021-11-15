#ifndef _F4D_INSTANCEPARSER_HPP_
#define _F4D_INSTANCEPARSER_HPP_

#include <string.h>

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
    class InstanceParser {

        Flat2DVector<float> costs_matrix;
        std::vector<std::vector<int>> neighbors;

    protected:
        const std::string& path;
        int vehicle_capacity;
        std::vector<float> x_coordinates;
        std::vector<float> y_coordinates;
        std::vector<int> demands;

    public:
        explicit InstanceParser(const std::string& path_) : path(path_) { }

        bool parse() {

            FILE* file = fopen(path.c_str(), "r");

            if (!file) {
                return false;
            }

            char name[64];
            char edgeWeightType[64];
            char aux[200];

            if (fscanf(file, "NAME : %s\n", name) != 1) return false;
            if (fscanf(file, "COMMENT : %*[^\n]\n") != 0) return false;
            if (fscanf(file, "TYPE : %*[^\n]\n") != 0) return false;

            int matrix_size;
            if (fscanf(file, "DIMENSION : %d\n", &matrix_size) != 1) return false;
            if (fscanf(file, "EDGE_WEIGHT_TYPE : %s\n", edgeWeightType) != 1) return false;

            bool is_explicit = strcmp(edgeWeightType, "EUC_2D");

            if (is_explicit) {
                if (fscanf(file, "EDGE_WEIGHT_FORMAT : %s\n", aux) != 1) return false;
                if (fscanf(file, "NODE_COORD_TYPE : %*[^\n]\n") != 0) return false;
            }

            if (fscanf(file, "CAPACITY : %d\n", &vehicle_capacity) != 1) return false;

            costs_matrix.resize(matrix_size, matrix_size);

            this->x_coordinates.resize(matrix_size);
            this->y_coordinates.resize(matrix_size);

            if (!is_explicit) {
                if (fscanf(file, "NODE_COORD_SECTION\n") != 0) return false;

                int vertex_index;
                for (int i = 0; i < matrix_size; ++i) {
                    if (fscanf(file, "%d %f %f", &vertex_index, &x_coordinates[i], &y_coordinates[i]) != 3) return false;
                }

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

            } else {
                if (fscanf(file, "EDGE_WEIGHT_SECTION\n") != 0) return false;

                for (int i = 0; i < matrix_size; i++) {
                    for (int j = 0; j < matrix_size; j++) {
                        if (i == j) break;

                        float c_ij;
                        if (fscanf(file, "%f", &c_ij) != 1) return false;
                        costs_matrix.at(i, j, c_ij);
                        costs_matrix.at(j, i, costs_matrix.at(i, j));
                    }
                }

                if (fscanf(file, "\n") != 0) return false;
                if (fscanf(file, "NODE_COORD_SECTION\n") != 0) return false;
                int vertex_index;
                for (int i = 0; i < matrix_size; ++i) {
                    if (fscanf(file, "%d %f %f", &vertex_index, &x_coordinates[i], &y_coordinates[i]) != 3) return false;
                }
            }

            if (fscanf(file, "\n") != 0) return false;
            if (fscanf(file, "DEMAND_SECTION\n") != 0) return false;

            int vertex_index;
            demands.resize(matrix_size);
            for (int i = 0; i < matrix_size; i++) {
                if (fscanf(file, "%d %d", &vertex_index, &demands[i]) != 2) return false;
            }

            fclose(file);

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

            return true;
        }

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


}  // namespace cobra


#endif