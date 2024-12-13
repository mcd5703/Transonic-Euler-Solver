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
// Include the solver header
#include "solver.hpp"

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

    int NACA = config::DEFAULT_NACA;
    double Minf = config::DEFAULT_MINF;
    double alpha = config::DEFAULT_ALPHA;
    double gamInf = config::DEFAULT_GAMINF;
    int itmax = config::DEFAULT_ITMAX;
    double CFL = config::DEFAULT_CFL;
    int nDisp = config::DEFAULT_NDISP;
    int NJ = config::DEFAULT_GRID_NJ;
    int NK = config::DEFAULT_GRID_NK;
    double R2 = config::DEFAULT_GRID_R2;


    Solver(NACA, Minf, alpha, gamInf, itmax, CFL, nDisp, NJ, NK, R2);

	return 0;
}