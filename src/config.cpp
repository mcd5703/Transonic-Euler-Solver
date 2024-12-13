#include "config.h"
#include <iostream>
#include <fstream>
#include <map>

// Default input values
int DEFAULT_NACA        = 2412; // NACA Configuraton
double DEFAULT_MINF     = 0.82; // Freestream Mach number
double DEFAULT_ALPHA    = 0.0;  // Angle of attack
double DEFAULT_GAMINF   = 1.4;  // Isentropic Compression Ratio
int DEFAULT_ITMAX       = 3000; // Max iterations
double DEFAULT_CFL      = 50.0; // CFL Number
int DEFAULT_NDISP       = 25;   // Dispute every __ frames
int DEFAULT_GRID_NJ     = 41;   // Grid NJ
int DEFAULT_GRID_NK     = 81;   // Grid NK
int DEFAULT_GRID_R2     = 100;  // Grid R2

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

    DEFAULT_NACA = config.count("DEFAULT_NACA") ? std::stoi(config.at("DEFAULT_NACA")) : DEFAULT_NACA;
    DEFAULT_MINF = config.count("DEFAULT_MINF") ? std::stoi(config.at("DEFAULT_MINF")) : DEFAULT_MINF;
    DEFAULT_ALPHA = config.count("DEFAULT_ALPHA") ? std::stoi(config.at("DEFAULT_ALPHA")) : DEFAULT_ALPHA;
    DEFAULT_GAMINF = config.count("DEFAULT_GAMINF") ? std::stoi(config.at("DEFAULT_GAMINF")) : DEFAULT_GAMINF;
    DEFAULT_ITMAX = config.count("DEFAULT_ITMAX") ? std::stoi(config.at("DEFAULT_ITMAX")) : DEFAULT_ITMAX;
    DEFAULT_CFL = config.count("DEFAULT_CFL") ? std::stoi(config.at("DEFAULT_CFL")) : DEFAULT_CFL;
    DEFAULT_NDISP = config.count("DEFAULT_NDISP") ? std::stoi(config.at("DEFAULT_NDISP")) : DEFAULT_NDISP;

    DEFAULT_GRID_NJ = config.count("DEFAULT_GRID_NJ") ? std::stoi(config.at("DEFAULT_GRID_NJ")) : DEFAULT_GRID_NJ;
    DEFAULT_GRID_NK = config.count("DEFAULT_GRID_NK") ? std::stoi(config.at("DEFAULT_GRID_NK")) : DEFAULT_GRID_NK;
    DEFAULT_GRID_R2 = config.count("DEFAULT_GRID_R2") ? std::stoi(config.at("DEFAULT_GRID_R2")) : DEFAULT_GRID_R2;
};
