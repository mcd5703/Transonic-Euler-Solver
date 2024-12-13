#ifndef QDOTQ_HPP
#define QDOTQ_HPP

#include <vector>

inline double qdotq(
    const std::vector<std::vector<std::vector<double>>> &q1,
    const std::vector<std::vector<std::vector<double>>> &q2,
    int NJ, int NK
) {
    double H = 0.0;
    for (int j = 1; j < NJ - 1; j++) {
        for (int k = 0; k < NK; k++) {
            H += q1[j][k][0] * q2[j][k][0];
        }
    }
    return H;
}

#endif // QDOTQ_HPP
