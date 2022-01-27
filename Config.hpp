#ifndef _F4D_CONFIG_HPP_
#define _F4D_CONFIG_HPP_

#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

class Config {

public:
    Config(const std::string& path) {

        std::ifstream file(path);

        if (file.is_open()) {
            std::string line;
            while (std::getline(file, line)) {
                if (line.empty() || line.front() == '#') continue;
                auto [key, value] = split(line);
                data[key] = value;
            }
            file.close();
        }
    }

    bool get_bool(const std::string& key) {
        return get_int(key);
    }

    int get_int(const std::string& key) {
        return std::stoi(data[key]);
    }

    std::string get_string(const std::string& key) {
        return data[key];
    }

    float get_real(const std::string& key) {
        return std::stof(data[key]);
    }

private:
    std::unordered_map<std::string, std::string> data = {
        
        {"small_cache", "50"},
        {"small_tier2", "0"},
        {"small_fastopt", "0"},
        {"small_routemin", "1000"},
        {"small_granular_neighbors", "25"},
        
        {"medium_cache", "50"},
        {"medium_tier2", "0"},
        {"medium_fastopt", "0"},
        {"medium_routemin", "1000"},
        {"medium_granular_neighbors", "25"},

        {"large_cache", "50"},
        {"large_tier2", "0"},
        {"large_fastopt", "0"},
        {"large_routemin", "1000"},
        {"large_granular_neighbors", "25"},

        {"xlarge_cache", "50"},
        {"xlarge_tier2", "0"},
        {"xlarge_fastopt", "0"},
        {"xlarge_routemin", "1000"},
        {"xlarge_granular_neighbors", "25"},

        {"xxlarge_cache", "50"},
        {"xxlarge_tier2", "0"},
        {"xxlarge_fastopt", "0"},
        {"xxlarge_routemin", "1000"},
        {"xxlarge_granular_neighbors", "25"},

    };


    std::pair<std::string, std::string> split(const std::string& line) {

        std::vector<std::string> tokens;

        size_t i = 0;
        std::string key;
        while (line[i] != '=') key += line[i++];

        ++i;  // skip =

        std::string value;
        while (i < line.size()) value += line[i++];

        return {key, value};
    }
};

#endif