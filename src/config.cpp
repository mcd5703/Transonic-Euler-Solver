#include "config.h"
#include <iostream>
#include <fstream>
#include <map>

int config::TIME_GRABBING_TOOL = 100;
int config::TIME_LANDING_PROCESS = 100;
int config::TIME_USING_TOOL = 100;
std::string config::DATA_FILE = "data.csv";

void config::parseConfigFile(std::string fqp)
{
    std::map<std::string, std::string> config;
    std::ifstream file(fqp);
    std::string line, key, value;
    while (std::getline(file, line))
    {
        std::size_t pos = line.find(' ');
        if (std::string::npos != pos)
        {
            key = line.substr(0, pos);
            value = line.substr(pos + 1, (line.size() - pos) - 2);
            config[key] = value;
        }
    }
    for (auto const &[key, value] : config)
        std::cout << key << "\t" << value << std::endl;

    TIME_LANDING_PROCESS = config.count("TIME_LANDING_PROCESS") ? std::stoi(config.at("TIME_LANDING_PROCESS")) : TIME_LANDING_PROCESS;
    TIME_GRABBING_TOOL = config.count("TIME_GRABBING_TOOL") ? std::stoi(config.at("TIME_GRABBING_TOOL")) : TIME_GRABBING_TOOL;
    TIME_USING_TOOL = config.count("TIME_USING_TOOL") ? std::stoi(config.at("TIME_USING_TOOL")) : TIME_USING_TOOL;
    DATA_FILE = config.count("DATA_FILE") ? config.at("DATA_FILE") : DATA_FILE;
};
