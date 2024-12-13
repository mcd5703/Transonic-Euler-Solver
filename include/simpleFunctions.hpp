#ifndef SIMPLEFUNCTIONS_HPP
#define SIMPLEFUNCTIONS_HPP

#include <vector>

// qdotq computes the inner product of q1 and q2 arrays over specified indices.
double qdotq(
    const std::vector<std::vector<std::vector<double>>> &q1,
    const std::vector<std::vector<std::vector<double>>> &q2,
    int NJ, int NK
);

// fluxjac computes the product A*v1, where A is derived from q, xx, xy, gamInf.
void fluxjac(
    const std::vector<double> &v1,
    const std::vector<double> &q,
    double xx, double xy, double gamInf,
    std::vector<double> &fout
);

#endif // SIMPLEFUNCTIONS_HPP
