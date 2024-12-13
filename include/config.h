#pragma once
#include <string>

namespace config
{
    extern int DEFAULT_NACA;
    extern double DEFAULT_MINF;
    extern double DEFAULT_ALPHA;
    extern double DEFAULT_GAMINF;
    extern int DEFAULT_ITMAX;
    extern int DEFAULT_CFL;
    extern int DEFAULT_NDISP;
    extern int DEFAULT_GRID_NJ;
    extern int DEFAULT_GRID_NK;
    extern int DEFAULT_GRID_R2;

    const static inline std::string RESULTS_DIR = "results/";

    /// @brief Quick (i.e., dumb) configuration file parser so we only need to compile the
    /// code one time and execute with different parameters
    /// @param fqp path to the configuration file
    void parseConfigFile(std::string fqp);
}