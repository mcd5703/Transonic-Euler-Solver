#include <Eigen/Dense>
#include <vector>

// Use Eigen for matrices
using namespace Eigen;
typedef MatrixXd Matrix;

// Function: metrics
void metrics(const Matrix &x, const Matrix &y, Matrix &xx, Matrix &xy, Matrix &yx, Matrix &yy, Matrix &vol, int NJ, int NK) {
    Matrix xxsi = Matrix::Zero(NJ, NK);
    Matrix yxsi = Matrix::Zero(NJ, NK);
    Matrix xeta = Matrix::Zero(NJ, NK);
    Matrix yeta = Matrix::Zero(NJ, NK);

    // Calculate xxsi and yxsi (central differences in J direction)
    for (int k = 0; k < NK; ++k) {
        xxsi(0, k) = x(1, k) - x(0, k);
        yxsi(0, k) = y(1, k) - y(0, k);
        for (int j = 1; j < NJ - 1; ++j) {
            xxsi(j, k) = 0.5 * (x(j + 1, k) - x(j - 1, k));
            yxsi(j, k) = 0.5 * (y(j + 1, k) - y(j - 1, k));
        }
        xxsi(NJ - 1, k) = x(NJ - 1, k) - x(NJ - 2, k);
        yxsi(NJ - 1, k) = y(NJ - 1, k) - y(NJ - 2, k);
    }

    // Calculate xeta and yeta (central differences in K direction with periodicity)
    for (int j = 0; j < NJ; ++j) {
        xeta(j, 0) = 0.5 * (x(j, 1) - x(j, NK - 2)); // Periodic boundary
        yeta(j, 0) = 0.5 * (y(j, 1) - y(j, NK - 2));
        for (int k = 1; k < NK - 1; ++k) {
            xeta(j, k) = 0.5 * (x(j, k + 1) - x(j, k - 1));
            yeta(j, k) = 0.5 * (y(j, k + 1) - y(j, k - 1));
        }
        xeta(j, NK - 1) = xeta(j, 0); // Periodicity
        yeta(j, NK - 1) = yeta(j, 0);
    }

    // Calculate grid metrics: vol, xx, xy, yx, yy
    for (int j = 0; j < NJ; ++j) {
        for (int k = 0; k < NK; ++k) {
            vol(j, k) = xxsi(j, k) * yeta(j, k) - xeta(j, k) * yxsi(j, k);
            xx(j, k) = yeta(j, k);
            xy(j, k) = -xeta(j, k);
            yx(j, k) = -yxsi(j, k);
            yy(j, k) = xxsi(j, k);
        }
    }
}
