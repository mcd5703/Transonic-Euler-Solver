#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined _MSC_VER
#include <direct.h>
#elif defined __GNUC__
#include <sys/types.h>
#include <sys/stat.h>
#endif

#define _CRT_INTERNAL_NONSTDC_NAMES 1
#include <sys/stat.h>
#if !defined(S_ISREG) && defined(S_IFMT) && defined(S_IFREG)
#define S_ISREG(m) (((m) & S_IFMT) == S_IFREG)
#endif
#if !defined(S_ISDIR) && defined(S_IFMT) && defined(S_IFDIR)
#define S_ISDIR(m) (((m) & S_IFMT) == S_IFDIR)
#endif

/// @brief A quick-and-dirty method to make a directory safely
/// @param dir_name path to the directory to be made
void create_directory_if_not_exists(const char *dir_name)
{
    struct stat st;
    // Check if the directory exists
    if (stat(dir_name, &st) == -1)
    {
// Directory does not exist, create it
#if defined _MSC_VER
        if (_mkdir(dir_name) == 0)
#elif defined __GNUC__
        if (mkdir(dir_name, 0755) == 0)
#endif

            printf("directory '%s' created successfully.\n", dir_name);
        else
            perror("error creating directory");
    }
    else if (S_ISDIR(st.st_mode))
        printf("directory '%s' already exists.\n", dir_name);
    else
        printf("a file with the name '%s' exists, but it is not a directory.\n", dir_name);
}