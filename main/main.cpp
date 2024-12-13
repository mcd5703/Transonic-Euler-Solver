/** *******************************************************************
 * @file
 * 
 * @brief main source file
 ******************************************************************** */
#include <iostream>
// #include "Question1.h"
#include "config.h"
#ifdef _WIN32
	#include <windows.h>
	#include "getopt.h"
#else
	#include <unistd.h>
#endif

// Provided by Dr. Simon Miller
void parseArgs(int argc, char **argv)
{
	
	std::string fqp;
	// Loop through all user-provided arguments and process the options
	char c;
	while ((c = getopt(argc, argv, "c:h")) != -1)
	{
		switch (c)
		{
		case 'c':
			fqp = optarg;
			break;
		case 'h':
		default:
			printf("main [options] \n");
			printf("-c    Configuration File Path\n");
			exit(-1);
		}
	}
	if (!fqp.empty())
		config::parseConfigFile(fqp);
}

int main(int argc, char *argv[])
{
	parseArgs(argc, argv);

    // Output the values for each variable in "config.h"
	std::cout << "-----------------------------------------------------------" << std::endl;
	std::cout << "NOTICE: The following defaults were read through config.inp... " << std::endl;
    std::cout << "DEFAULT_NACA: " << config::DEFAULT_NACA << std::endl;
    std::cout << "DEFAULT_MINF: " << config::DEFAULT_MINF << std::endl;
    std::cout << "DEFAULT_ALPHA: " << config::DEFAULT_ALPHA << std::endl;
    std::cout << "DEFAULT_GAMINF: " << config::DEFAULT_GAMINF << std::endl;
    std::cout << "DEFAULT_ITMAX: " << config::DEFAULT_ITMAX << std::endl;
    std::cout << "DEFAULT_CFL: " << config::DEFAULT_CFL << std::endl;
    std::cout << "DEFAULT_NDISP: " << config::DEFAULT_NDISP << std::endl;
    std::cout << "DEFAULT_GRID_NJ: " << config::DEFAULT_GRID_NJ << std::endl;
    std::cout << "DEFAULT_GRID_NK: " << config::DEFAULT_GRID_NK << std::endl;
    std::cout << "DEFAULT_GRID_R2: " << config::DEFAULT_GRID_R2 << std::endl;
	std::cout << "-----------------------------------------------------------" << std::endl;

	return 0;
}