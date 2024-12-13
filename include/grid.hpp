#ifndef GRID_HPP
#define GRID_HPP

#include <vector>

// Prototype for makegrid function
// Parameters:
//  NACA, R2, NJ, NK: integers/doubles to define the geometry and grid resolution
//  x, y: 2D vectors to be filled with the generated grid coordinates
void makegrid(int NACA, double R2, int NJ, int NK,
                std::vector<std::vector<double>> &x,
                std::vector<std::vector<double>> &y);

// Prototype for metrics function
// Parameters:
//  x, y: the coordinate arrays from makegrid
//  NJ, NK: grid dimensions
//  xx, xy, yx, yy, vol: metrics and volume arrays to be computed
void metrics(const std::vector<std::vector<double>> &x,
                const std::vector<std::vector<double>> &y,
                int NJ, int NK,
                std::vector<std::vector<double>> &xx,
                std::vector<std::vector<double>> &xy,
                std::vector<std::vector<double>> &yx,
                std::vector<std::vector<double>> &yy,
                std::vector<std::vector<double>> &vol);

#endif // GRID_HPP
