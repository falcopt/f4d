#include "arg_parser.hpp"

void print_help() {

    std::cout << "Usage: filo <path-to-instance> [OPTIONS]\n\n";

    std::cout << "Available options\n";

    std::cout << TOKEN_OUTPATH << " STRING\t\tOutput directory (default: " << DEFAULT_OUTPATH << ")\n";
    std::cout << TOKEN_PARSER << " STRING\t\t\tParser type, it can be X, K or Z (default: " << DEFAULT_PARSER << ")\n";
    std::cout << TOKEN_TOLERANCE << " FLOAT\t\tFloating point tolerance (default: " << DEFAULT_TOLERANCE << ")\n";
    std::cout << TOKEN_SPARSIFICATION_RULE1_NEIGHBORS
              << " INT\tNeighbors per vertex in granular neighborhoods (default: " << DEFAULT_SPARSIFICATION_RULE1_NEIGHBORS << ")\n";
    std::cout << TOKEN_SOLUTION_CACHE_HISTORY << " INT\t\t\tSelective cache dimension (default: " << DEFAULT_SOLUTION_CACHE_HISTORY << ")\n";
    std::cout << TOKEN_ROUTEMIN_ITERATIONS << " INT\tMax route minimization iterations (default: " << DEFAULT_ROUTEMIN_ITERATIONS << ")\n";
    std::cout << TOKEN_COREOPT_ITERATIONS << " INT\tCore optimization iterations (default: " << DEFAULT_COREOPT_ITERATIONS << ")\n";
    std::cout << TOKEN_SPARSIFICATION_FACTOR << " FLOAT\tInitial sparsification factor gamma base (default: " << DEFAULT_SPARSIFICATION_FACTOR << ")\n";
    std::cout << TOKEN_SPARSIFICATION_MULTIPLIER << " FLOAT\t\tGranular reduction factor delta (default: " << DEFAULT_SPARSIFICATION_MULTIPLIER << ")\n";
    std::cout << TOKEN_SHAKING_LB_FACTOR << " FLOAT\tShaking lower bound factor (default: " << DEFAULT_SHAKING_LB_FACTOR << ")\n";
    std::cout << TOKEN_SHAKING_UB_FACTOR << " FLOAT\tShaking upper bound factor (default: " << DEFAULT_SHAKING_UB_FACTOR << ")\n";
    std::cout << TOKEN_SEED << " INT\t\t\tSeed (default: " << DEFAULT_SEED << ")\n";
    std::cout << TOKEN_ROUND << " 1-0\t\t\tRound costs (default: " << DEFAULT_ROUND << ")\n";

    std::cout << "\nReport bugs to luca.accorsi4 or f.cavaliere and the domain is unibo.it\n";
}

Parameters parse_command_line_arguments(int argc, char* argv[]) {

    if (argc == 1) {
        std::cout << "Error: missing input instance.\n\n";
        print_help();
        exit(EXIT_FAILURE);
    }

    if (argc == 2 && std::string(argv[1]) == TOKEN_HELP) {
        print_help();
        exit(EXIT_SUCCESS);
    }

    auto parameters = Parameters(std::string(argv[1]));

    // dimacs run
    if (argc == 3) {

        parameters.set(TOKEN_ROUND, argv[2]);
        return parameters;
    }

    for (auto n = 2; n < argc; n += 2) {

        auto token = std::string(argv[n]);

        if (token == TOKEN_HELP) {
            print_help();
            exit(EXIT_SUCCESS);
        }

        if (n + 1 >= argc) {
            std::cout << "Error: missing value for '" << token << "'.\n\n";
            print_help();
            exit(EXIT_FAILURE);
        }
        auto value = std::string(argv[n + 1]);

        parameters.set(token, value);
    }

    return parameters;
}