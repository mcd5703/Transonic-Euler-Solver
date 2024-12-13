#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>

void Solver(int NACA, double Minf, double alpha, double gamInf,
            int itmax, double CFL, int nDisp, int NJ, int NK, double R2);

// Prototype for the visualizeGrid function
// This will just plot the x,y grid using matplot++
void visualizeGrid(const std::vector<std::vector<double>> &x,
                    const std::vector<std::vector<double>> &y,
                    int NJ, int NK);

#endif // SOLVER_HPP
