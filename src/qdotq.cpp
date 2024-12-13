#include "simpleFunctions.hpp"

double qdotq(
    const std::vector<std::vector<std::vector<double>>> &q1,
    const std::vector<std::vector<std::vector<double>>> &q2,
    int NJ, int NK
)
{
    double H = 0.0;
    // j=1 to NJ-2 in C++ (since MATLAB's 2 to NJ-1 is shifted by one for 0-based)
    for (int j = 1; j < NJ-1; j++) {
        for (int k = 0; k < NK; k++) {
            H += q1[j][k][0]*q2[j][k][0];
        }
    }
    return H;
}
