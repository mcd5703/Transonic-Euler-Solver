#include <Eigen/Dense>
#include <vector>

// Use Eigen for matrices
using namespace Eigen;
typedef MatrixXd Matrix;

// Function: qdotq
double qdotq(const std::vector<Matrix> &q1, const std::vector<Matrix> &q2, int NJ, int NK) {
    double H = 0.0;

    // Sum over the valid ranges
    for (int j = 1; j < NJ - 1; ++j) { // MATLAB 2:NJ-1 -> C++ 1:NJ-2
        for (int k = 0; k < NK; ++k) { // MATLAB 1:NK -> C++ 0:NK-1
            H += q1[0](j, k) * q2[0](j, k); // First component of q1 and q2
        }
    }

    return H;
}
