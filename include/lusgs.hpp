#ifndef LUSGS_HPP
#define LUSGS_HPP

#include <vector>

// lusgs function performs a LU-SGS relaxation step given R, q, metrics, dt, etc.
// Returns dq array after the procedure.
void lusgs(
    const std::vector<std::vector<std::vector<double>>> &R,
    const std::vector<std::vector<std::vector<double>>> &q,
    const std::vector<std::vector<double>> &xx,
    const std::vector<std::vector<double>> &xy,
    const std::vector<std::vector<double>> &yx,
    const std::vector<std::vector<double>> &yy,
    const std::vector<std::vector<double>> &vol,
    const std::vector<std::vector<double>> &dt,
    double gamInf, int NJ, int NK,
    std::vector<std::vector<std::vector<double>>> &dq
);

#endif // LUSGS_HPP
