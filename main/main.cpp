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


	return 0;
}