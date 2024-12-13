#ifndef RESID_HPP
#define RESID_HPP

#include <vector>

// resid function computes the residual R and time step dt given the current solution q and other parameters.
// q, xx, xy, yx, yy, vol, qInfPrim are input.
// R and dt are output parameters passed by reference.
void resid(
    const std::vector<std::vector<std::vector<double>>> &q,
    const std::vector<std::vector<double>> &xx,
    const std::vector<std::vector<double>> &xy,
    const std::vector<std::vector<double>> &yx,
    const std::vector<std::vector<double>> &yy,
    const std::vector<std::vector<double>> &vol,
    const std::vector<double> &qInfPrim,
    int NJ, int NK, double CFL,
    std::vector<std::vector<std::vector<double>>> &R,
    std::vector<std::vector<double>> &dt
);

#endif // RESID_HPP
