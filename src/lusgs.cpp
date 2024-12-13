#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <algorithm>

// Use Eigen for matrices
using namespace Eigen;
typedef MatrixXd Matrix;

// Function prototype for fluxjac (to be implemented separately)
Vector4d fluxjac(const Vector4d &dq, const Vector4d &q, double xx_vol, double xy_vol, double gamInf);

// LU-SGS function
std::vector<Matrix> lusgs(const std::vector<Matrix> &R, const std::vector<Matrix> &q, 
                        const Matrix &xx, const Matrix &xy, const Matrix &yx, const Matrix &yy, 
                        const Matrix &vol, const Matrix &dt, double gamInf, int NJ, int NK) {
    // Initialize dq and dqs as zero 3D arrays
    std::vector<Matrix> dq(4, Matrix::Zero(NJ, NK));
    std::vector<Matrix> dqs(4, Matrix::Zero(NJ, NK));

    Matrix sigj = Matrix::Zero(NJ, NK);
    Matrix sigk = Matrix::Zero(NJ, NK);

    // Build spectral radii
    for (int j = 0; j < NJ; ++j) {
        for (int k = 0; k < NK; ++k) {
            double rho = q[0](j, k);
            double u = q[1](j, k) / rho;
            double v = q[2](j, k) / rho;
            double T = (gamInf - 1.0) * (q[3](j, k) / rho - 0.5 * (u * u + v * v));
            sigj(j, k) = (fabs(xx(j, k) * u + xy(j, k) * v) + sqrt((xx(j, k) * xx(j, k) + xy(j, k) * xy(j, k)) * gamInf * T)) / vol(j, k);
            sigk(j, k) = (fabs(yx(j, k) * u + yy(j, k) * v) + sqrt((yx(j, k) * yx(j, k) + yy(j, k) * yy(j, k)) * gamInf * T)) / vol(j, k);
        }
    }

    Vector4d tmp1, tmp2;

    // Forward substitution loop
    for (int j = 0; j < NJ; ++j) {
        for (int k = 0; k < NK; ++k) {
            int jm = std::max(j - 1, 0);
            int km = std::max(k - 1, 0);

            double facx = (j > 0) ? 1.0 : 2.0;
            double facy = (k > 0) ? 1.0 : 2.0;

            tmp1.setZero();
            tmp2.setZero();

            if (j > 0) {
                tmp1 = fluxjac(dqs[jm][k], q[jm][k], xx(jm, k) / vol(jm, k), xy(jm, k) / vol(jm, k), gamInf);
                tmp1 = 0.5 * (tmp1 + sigj(jm, k) * dqs[jm][k]);
            }

            if (k > 0) {
                tmp2 = fluxjac(dqs[j][km], q[j][km], yx(j, km) / vol(j, km), yy(j, km) / vol(j, km), gamInf);
                tmp2 = 0.5 * (tmp2 + sigk(j, km) * dqs[j][km]);
            }

            double D = 1.0 + dt(j, k) * (facx * sigj(j, k) + facy * sigk(j, k));

            for (int n = 0; n < 4; ++n) {
                dqs[n](j, k) = (-R[n](j, k) * D * vol(j, k) + dt(j, k) * (tmp1(n) + tmp2(n))) / D;
            }
        }
    }

    // Clean up periodicity
    for (int j = 0; j < NJ; ++j) {
        for (int n = 0; n < 4; ++n) {
            dqs[n](j, 0) = 0.5 * (dqs[n](j, 0) + dqs[n](j, NK - 1));
            dqs[n](j, NK - 1) = dqs[n](j, 0);
        }
    }

    // Backward substitution loop
    for (int j = NJ - 1; j >= 0; --j) {
        for (int k = NK - 1; k >= 0; --k) {
            int jp = std::min(j + 1, NJ - 1);
            int kp = std::min(k + 1, NK - 1);

            double facx = (j < NJ - 1) ? 1.0 : 2.0;
            double facy = (k < NK - 1) ? 1.0 : 2.0;

            tmp1.setZero();
            tmp2.setZero();

            if (j < NJ - 1) {
                tmp1 = fluxjac(dq[jp][k], q[jp][k], xx(jp, k) / vol(jp, k), xy(jp, k) / vol(jp, k), gamInf);
                tmp1 = 0.5 * (tmp1 - sigj(jp, k) * dq[jp][k]);
            }

            if (k < NK - 1) {
                tmp2 = fluxjac(dq[j][kp], q[j][kp], yx(j, kp) / vol(j, kp), yy(j, kp) / vol(j, kp), gamInf);
                tmp2 = 0.5 * (tmp2 - sigk(j, kp) * dq[j][kp]);
            }

            double D = 1.0 + dt(j, k) * (facx * sigj(j, k) + facy * sigk(j, k));

            for (int n = 0; n < 4; ++n) {
                dq[n](j, k) = (dqs[n](j, k) - dt(j, k) * (tmp1(n) + tmp2(n))) / D;
            }
        }
    }

    // Clean up periodicity
    for (int j = 0; j < NJ; ++j) {
        for (int n = 0; n < 4; ++n) {
            dq[n](j, 0) = 0.5 * (dq[n](j, 0) + dq[n](j, NK - 1));
            dq[n](j, NK - 1) = dq[n](j, 0);
        }
    }

    // Divide by volume
    for (int n = 0; n < 4; ++n) {
        for (int j = 0; j < NJ; ++j) {
            for (int k = 0; k < NK; ++k) {
                dq[n](j, k) /= vol(j, k);
            }
        }
    }

    return dq;
}