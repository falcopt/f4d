#ifndef _F4D_ARGPARSER_HPP_
#define _F4D_ARGPARSER_HPP_

#include <filesystem>
#include <iostream>
#include <optional>
#include <string>
#include <utility>

/* Default parameters */
#define DEFAULT_OUTPATH ("./")
#define DEFAULT_PARSER ("TSPLIB95")
#define DEFAULT_SOLUTION_CACHE_HISTORY (50)
#define DEFAULT_CW_LAMBDA (1.0f)
#define DEFAULT_CW_NEIGHBORS (100)
#define DEFAULT_ROUTEMIN_ITERATIONS (1000)
#define DEFAULT_COREOPT_ITERATIONS (100000)
#define DEFAULT_FASTOPT_ITERATIONS (0)
#define DEFAULT_SPARSIFICATION_RULE1_NEIGHBORS (25)
#define DEFAULT_SPARSIFICATION_FACTOR (0.25f)
#define DEFAULT_SPARSIFICATION_MULTIPLIER (0.50f)
#define DEFAULT_SHAKING_LB_FACTOR (0.375f)
#define DEFAULT_SHAKING_UB_FACTOR (0.85f)
#define DEFAULT_TOLERANCE (0.01f)
#define DEFAULT_SEED (0)
#define DEFAULT_ROUND (1)
#define DEFAULT_BKS_VALUE (std::nullopt)
#define DEFAULT_SA_INIT_FACTOR (0.1f)
#define DEFAULT_SA_FINAL_FACTOR (0.001f)
#define DEFAULT_TIER2 (0)

/* Tokens */
#define TOKEN_OUTPATH ("--outpath")
#define TOKEN_PARSER ("--parser")
#define TOKEN_TOLERANCE ("--tolerance")
#define TOKEN_SPARSIFICATION_RULE1_NEIGHBORS ("--granular-neighbors")
#define TOKEN_SOLUTION_CACHE_HISTORY ("--cache")
#define TOKEN_ROUTEMIN_ITERATIONS ("--routemin-iterations")
#define TOKEN_COREOPT_ITERATIONS ("--coreopt-iterations")
#define TOKEN_FASTOPT_ITERATIONS ("--fastopt-iterations")
#define TOKEN_SPARSIFICATION_FACTOR ("--granular-gamma-base")
#define TOKEN_SPARSIFICATION_MULTIPLIER ("--granular-delta")
#define TOKEN_SHAKING_LB_FACTOR ("--shaking-lower-bound")
#define TOKEN_SHAKING_UB_FACTOR ("--shaking-upper-bound")
#define TOKEN_SEED ("--seed")
#define TOKEN_HELP ("--help")
#define TOKEN_ROUND ("--round")
#define TOKEN_BKS_VALUE ("--bks-value")
#define TOKEN_SA_INIT_FACTOR ("--sa-initial-factor")
#define TOKEN_SA_FINAL_FACTOR ("--sa-final-factor")
#define TOKEN_TIER2 ("--tier2")

class Parameters {

private:
    std::string instance_path;
    std::string outpath = DEFAULT_OUTPATH;
    std::string parser = DEFAULT_PARSER;
    float tolerance = DEFAULT_TOLERANCE;
    int solution_cache_history = DEFAULT_SOLUTION_CACHE_HISTORY;
    float cw_lambda = DEFAULT_CW_LAMBDA;
    int cw_neighbors = DEFAULT_CW_NEIGHBORS;
    int routemin_iterations = DEFAULT_ROUTEMIN_ITERATIONS;
    int coreopt_iterations = DEFAULT_COREOPT_ITERATIONS;
    int fastopt_iterations = DEFAULT_FASTOPT_ITERATIONS;
    int sparsification_rule_neighbors = DEFAULT_SPARSIFICATION_RULE1_NEIGHBORS;
    float gamma_base = DEFAULT_SPARSIFICATION_FACTOR;
    float delta = DEFAULT_SPARSIFICATION_MULTIPLIER;
    float shaking_lb_factor = DEFAULT_SHAKING_LB_FACTOR;
    float shaking_ub_factor = DEFAULT_SHAKING_UB_FACTOR;
    int seed = DEFAULT_SEED;
    bool round_costs = DEFAULT_ROUND;
    std::optional<float> bks_value = DEFAULT_BKS_VALUE;
    float sa_initial_factor = DEFAULT_SA_INIT_FACTOR;
    float sa_final_factor = DEFAULT_SA_FINAL_FACTOR;
    int tier2 = DEFAULT_TIER2;

public:
    explicit Parameters(std::string instance_path_) : instance_path(std::move(instance_path_)) { }

    inline int get_solution_cache_size() const {
        return solution_cache_history;
    }
    inline float get_cw_lambda() const {
        return cw_lambda;
    }
    inline int get_cw_neighbors() const {
        return cw_neighbors;
    }
    inline int get_routemin_iterations() const {
        return routemin_iterations;
    }
    inline int get_coreopt_iterations() const {
        return coreopt_iterations;
    }
    inline int get_fastopt_iterations() const {
        return fastopt_iterations;
    }
    inline int get_sparsification_rule_neighbors() const {
        return sparsification_rule_neighbors;
    }
    inline float get_gamma_base() const {
        return gamma_base;
    }
    inline float get_delta() const {
        return delta;
    }
    inline float get_shaking_lb_factor() const {
        return shaking_lb_factor;
    }
    inline float get_shaking_ub_factor() const {
        return shaking_ub_factor;
    }
    inline float get_tolerance() const {
        return tolerance;
    }
    inline std::string get_instance_path() const {
        return instance_path;
    }
    inline std::string get_outpath() const {
        return outpath;
    }
    inline std::string get_parser() const {
        return parser;
    }
    inline int get_seed() const {
        return seed;
    }
    inline bool needs_round_costs() const {
        return round_costs;
    }
    inline auto get_bks_value() const {
        return bks_value;
    }
    inline auto get_sa_initial_factor() const {
        return sa_initial_factor;
    }
    inline auto get_sa_final_factor() const {
        return sa_final_factor;
    }

    inline auto get_tier2() const {
        return tier2;
    }

    void set(const std::string& key, const std::string& value) {

        if (key == TOKEN_OUTPATH) {
            outpath = value;
            if (outpath.back() != std::filesystem::path::preferred_separator) {
                outpath += std::filesystem::path::preferred_separator;
            }
        } else if (key == TOKEN_PARSER) {
            parser = value;
        } else if (key == TOKEN_TOLERANCE) {
            tolerance = std::stof(value);
        } else if (key == TOKEN_SPARSIFICATION_RULE1_NEIGHBORS) {
            sparsification_rule_neighbors = std::stoi(value);
        } else if (key == TOKEN_SOLUTION_CACHE_HISTORY) {
            solution_cache_history = std::stoi(value);
        } else if (key == TOKEN_ROUTEMIN_ITERATIONS) {
            routemin_iterations = std::stoi(value);
        } else if (key == TOKEN_COREOPT_ITERATIONS) {
            coreopt_iterations = std::stoi(value);
        } else if (key == TOKEN_FASTOPT_ITERATIONS) {
            fastopt_iterations = std::stoi(value);
        } else if (key == TOKEN_SPARSIFICATION_FACTOR) {
            gamma_base = std::stof(value);
        } else if (key == TOKEN_SPARSIFICATION_MULTIPLIER) {
            delta = std::stof(value);
        } else if (key == TOKEN_SHAKING_LB_FACTOR) {
            shaking_lb_factor = std::stof(value);
        } else if (key == TOKEN_SHAKING_UB_FACTOR) {
            shaking_ub_factor = std::stof(value);
        } else if (key == TOKEN_SEED) {
            seed = std::stoi(value);
        } else if (key == TOKEN_ROUND) {
            round_costs = std::stoi(value);
        } else if (key == TOKEN_BKS_VALUE) {
            bks_value = {std::stof(value)};
        } else if (key == TOKEN_SA_INIT_FACTOR) {
            sa_initial_factor = std::stof(value);
        } else if (key == TOKEN_SA_FINAL_FACTOR) {
            sa_final_factor = std::stof(value);
        } else if (key == TOKEN_TIER2) {
            tier2 = std::stoi(value);
        } else {
            std::cout << "Error: unknown argument '" << key << "'. Try --help for more information.\n";
            exit(EXIT_SUCCESS);
        }
    }
};

void print_help();
Parameters parse_command_line_arguments(char* argv[]);

#endif