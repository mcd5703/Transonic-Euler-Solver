#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>

// Use Eigen for matrices
using namespace Eigen;
typedef MatrixXd Matrix;

// Function to generate the grid
void makegrid(int NACA, double R2, int NJ, int NK, Matrix &x, Matrix &y) {
    x = Matrix::Zero(NJ, NK);
    y = Matrix::Zero(NJ, NK);

    // Extract NACA parameters
    int d1 = NACA / 1000;
    double m = d1 / 100.0;
    NACA -= 1000 * d1;

    int d2 = NACA / 100;
    double p = d2 / 10.0;
    NACA -= 100 * d2;

    double t = NACA / 100.0;
    double eps = 4 * sqrt(3) * t / 9;

    double z1 = 1.0;
    double z2 = z1 - 4 * (1 + 2 * eps) / (pow((1 + 2 * eps), 2) + 2 * (1 + 2 * eps) + 1);

    // Initialize variables
    std::vector<double> phi(NK), theta(NK), dtheta(NK), phiA(NK), error(NK);
    std::vector<std::complex<double>> zeta(NJ * NK), zeta0(NJ);

    for (int k = 0; k < NK; ++k) {
        phi[k] = 2.0 * M_PI * k / (NK - 1);
        theta[k] = phi[k];
    }

    double errorMax = 999.0;

    while (errorMax > 1e-10) {
        // Compute airfoil surface points
        for (int k = 0; k < NK; ++k) {
            x(0, k) = 0.5 * (1 + cos(theta[k]));
            y(0, k) = 5 * t * (0.2969 * sqrt(x(0, k)) - 0.1260 * x(0, k) - 0.3516 * x(0, k) * x(0, k) +
                               0.2843 * pow(x(0, k), 3) - 0.1036 * pow(x(0, k), 4));
            if (theta[k] > M_PI) y(0, k) = -y(0, k);

            if (x(0, k) < p) {
                y(0, k) += m * (2 * p * x(0, k) - x(0, k) * x(0, k)) / (p * p);
            } else {
                y(0, k) += m * ((1 - 2 * p) + 2 * p * x(0, k) - x(0, k) * x(0, k)) / pow((1 - p), 2);
            }
        }

        // Map the surface into the pseudo-circular plane
        std::complex<double> zeta1 = -0.5 * (z2 - z1), zeta2 = 0.5 * (z2 - z1);
        std::complex<double> sum = 0;

        for (int k = 0; k < NK; ++k) {
            std::complex<double> z = {x(0, k), y(0, k)};
            std::complex<double> RHS = pow((z - z1) / (z - z2), M_PI / (2 * M_PI - atan2(y(0, k), x(0, k))));
            zeta[k] = (zeta1 - RHS * zeta2) / (1.0 - RHS);

            sum += zeta[k];
        }

        std::complex<double> zeta0 = sum / (double)NK;

        // Compute phiA and errors
        for (int k = 0; k < NK; ++k) {
            phiA[k] = atan2(imag(zeta[k] - zeta0), real(zeta[k] - zeta0));
            if ((theta[k] - M_PI / 2 > 0) && (phiA[k] < 0)) phiA[k] += 2 * M_PI;
        }

        double phi0 = phiA[0];
        phiA[NK - 1] = phi0 + 2 * M_PI;

        for (int k = 0; k < NK; ++k) error[k] = phiA[k] - (phi[k] + phi0);
        errorMax = *std::max_element(error.begin(), error.end());

        // Update theta
        for (int k = 1; k < NK - 1; ++k) {
            double dphidtheta = (phiA[k + 1] - phiA[k - 1]) / (theta[k + 1] - theta[k - 1]);
            dtheta[k] = -error[k] / dphidtheta;
        }
        for (int k = 0; k < NK; ++k) theta[k] += dtheta[k];
    }

    // Generate pseudo-circular grid
    for (int k = 0; k < NK; ++k) {
        double R1 = abs(zeta[k] - zeta0);
        for (int j = 0; j < NJ; ++j) {
            double rad = R1 * R2 / (R2 + ((double)j / (NJ - 1)) * (R1 - R2));
            std::complex<double> z = zeta0 + rad * std::polar(1.0, phiA[k]);
            x(j, k) = real(z);
            y(j, k) = imag(z);
        }
    }
}
