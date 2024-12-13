#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <algorithm>

// Use Eigen for matrices and vectors
using namespace Eigen;
typedef MatrixXd Matrix;

// Function Prototypes
std::vector<double> HLLCFlux(const std::vector<double> &qL, const std::vector<double> &qR, 
                            double gamInf, double xA, double yA, double volA);

// Function: resid
std::pair<std::vector<Matrix>, Matrix> resid(const std::vector<Matrix> &q, 
                                            const Matrix &xx, const Matrix &xy, 
                                            const Matrix &yx, const Matrix &yy, 
                                            const Matrix &vol, const VectorXd &qInfPrim, 
                                            int NJ, int NK, double CFL) {
    // Residual and time step variables
    std::vector<Matrix> R(4, Matrix::Zero(NJ, NK));
    Matrix dt = Matrix::Zero(NJ, NK);

    // Free-stream parameters
    double rhoInf = qInfPrim(0);
    double uInf = qInfPrim(1);
    double vInf = qInfPrim(2);
    double TInf = qInfPrim(3);
    double gamInf = qInfPrim(4);
    double h0Inf = TInf * gamInf / (gamInf - 1.0) + 0.5 * (uInf * uInf + vInf * vInf);

    Matrix sigj = Matrix::Zero(NJ, NK), sigk = Matrix::Zero(NJ, NK);

    // Build spectral radii and CFL time step
    for (int j = 1; j < NJ - 1; ++j) {
        for (int k = 0; k < NK; ++k) {
            double rho = q[0](j, k);
            double u = q[1](j, k) / rho;
            double v = q[2](j, k) / rho;
            double T = (gamInf - 1.0) * (q[3](j, k) / rho - 0.5 * (u * u + v * v));
            sigj(j, k) = fabs(xx(j, k) * u + xy(j, k) * v) + sqrt((xx(j, k) * xx(j, k) + xy(j, k) * xy(j, k)) * gamInf * T);
            sigk(j, k) = fabs(yx(j, k) * u + yy(j, k) * v) + sqrt((yx(j, k) * yx(j, k) + yy(j, k) * yy(j, k)) * gamInf * T);

            double CFLeff = std::max(CFL / (1.0 + sqrt(vol(j, k))), 1.0);
            dt(j, k) = CFLeff * vol(j, k) / std::max(sigj(j, k), sigk(j, k));
        }
    }

    std::vector<double> qL(4), qR(4), fhat(4);

    // K residuals (periodic in K)
    for (int j = 1; j < NJ - 1; ++j) {
        for (int k = 0; k < NK - 1; ++k) {
            double yxA = 0.5 * (yx(j, k) + yx(j, k + 1));
            double yyA = 0.5 * (yy(j, k) + yy(j, k + 1));
            double volA = 0.5 * (vol(j, k) + vol(j, k + 1));

            // Periodic indices
            int km = (k == 0) ? NK - 1 : k - 1;
            int kp = k + 1;
            int kpp = (kp + 1 >= NK) ? 1 : kp + 1;

            // Left state reconstruction
            for (int n = 0; n < 4; ++n) {
                double r = (q[n](j, kp) - q[n](j, k)) / (q[n](j, k) - q[n](j, km) + 1e-8);
                double phi = std::max(0.0, std::min({2.0 * r, (1.0 + 2.0 * r) / 3.0, 2.0}));
                qL[n] = q[n](j, k) + 0.5 * phi * (q[n](j, k) - q[n](j, km));
            }

            // Right state reconstruction
            for (int n = 0; n < 4; ++n) {
                double r = (q[n](j, k) - q[n](j, kp)) / (q[n](j, kp) - q[n](j, kpp) + 1e-8);
                double phi = std::max(0.0, std::min({2.0 * r, (1.0 + 2.0 * r) / 3.0, 2.0}));
                qR[n] = q[n](j, kp) + 0.5 * phi * (q[n](j, kp) - q[n](j, kpp));
            }

            // Flux calculation
            fhat = HLLCFlux(qL, qR, gamInf, yxA, yyA, volA);

            // Subtract free-stream fluxes
            double UU = yxA * uInf + yyA * vInf;
            fhat[0] -= rhoInf * UU;
            fhat[1] -= rhoInf * UU * uInf + yxA * rhoInf * TInf;
            fhat[2] -= rhoInf * UU * vInf + yyA * rhoInf * TInf;
            fhat[3] -= rhoInf * UU * h0Inf;

            for (int n = 0; n < 4; ++n) {
                R[n](j, k) += fhat[n];
                R[n](j, kp) -= fhat[n];
            }
        }
    }

    // Volume-scale residuals and multiply by dt
    for (int j = 1; j < NJ - 1; ++j) {
        for (int k = 0; k < NK; ++k) {
            for (int n = 0; n < 4; ++n) {
                R[n](j, k) = dt(j, k) * R[n](j, k) / vol(j, k);
            }
        }
    }

    return {R, dt};
}

// Placeholder for HLLC flux calculation
std::vector<double> HLLCFlux(const std::vector<double> &qL, const std::vector<double> &qR, 
                            double gamInf, double xA, double yA, double volA) {
    std::vector<double> flux(4, 0.0);
    // Implement the HLLC flux solver here
    return flux;
}
