#pragma once
#include <string>

namespace config
{
    extern int TIME_LANDING_PROCESS;
    extern int TIME_GRABBING_TOOL;
    extern int TIME_USING_TOOL;
    extern std::string DATA_FILE;
    const static inline std::string RESULTS_DIR = "results/";

    /// @brief Quick (i.e., dumb) configuration file parser so we only need to compile the
    /// code one time and execute with different parameters
    /// @param fqp path to the configuration file
    void parseConfigFile(std::string fqp);
}