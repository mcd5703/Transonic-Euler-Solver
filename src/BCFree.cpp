#include <Eigen/Dense>
#include <vector>
#include <cmath>

// Use Eigen for matrices and vectors
using namespace Eigen;
typedef MatrixXd Matrix;
typedef std::vector<Matrix> Matrix3D; // 3D vector for qB (4 matrices for 4 state variables)

// Function: bcfree
Matrix3D bcfree(const Matrix3D &q, const VectorXd &qInfPrim, const Matrix &xx, const Matrix &xy, int NJ, int NK) {
    // Initialize output boundary condition matrix qB
    Matrix3D qB(4, Matrix::Zero(NJ, NK));

    // Free-stream quantities
    double rhoInf = qInfPrim(0);
    double uInf = qInfPrim(1);
    double vInf = qInfPrim(2);
    double TInf = qInfPrim(3);
    double gamInf = qInfPrim(4);

    double cInf = 1.0;  // Non-dimensionalized free-stream speed of sound

    // Loop over grid points
    for (int j = 0; j < NJ; ++j) {
        for (int k = 0; k < NK; ++k) {
            // Extract primitive variables
            double rho = q[0](j, k);
            double u = q[1](j, k) / rho;
            double v = q[2](j, k) / rho;
            double T = (gamInf - 1.0) * (q[3](j, k) / rho - 0.5 * (u * u + v * v));

            // Compute normal and tangential components of velocity
            double normFactor = std::sqrt(xx(j, k) * xx(j, k) + xy(j, k) * xy(j, k));
            double nx = xx(j, k) / normFactor;
            double ny = xy(j, k) / normFactor;

            double un = u * nx + v * ny;            // Normal velocity
            double unInf = uInf * nx + vInf * ny;  // Normal free-stream velocity
            double c = std::sqrt(gamInf * T);      // Local speed of sound

            double ut = u - un * nx;      // Tangential velocity in x
            double vt = v - un * ny;      // Tangential velocity in y
            double utInf = uInf - unInf * nx;
            double vtInf = vInf - unInf * ny;

            // Supersonic inflow: set to free stream
            if (unInf + cInf < 0) {
                qB[0](j, k) = rhoInf;
                qB[1](j, k) = rhoInf * uInf;
                qB[2](j, k) = rhoInf * vInf;
                qB[3](j, k) = rhoInf * (TInf / (gamInf - 1.0) + 0.5 * (uInf * uInf + vInf * vInf));
            }
            // Supersonic outflow: extrapolate
            else if (un > c) {
                for (int n = 0; n < 4; ++n) {
                    qB[n](j, k) = q[n](j, k);
                }
            }
            // Subsonic flow
            else {
                double R1 = un + 2.0 * c / (gamInf - 1.0);
                double R2 = unInf - 2.0 * cInf / (gamInf - 1.0);
                double s;

                if (un > 0) {  // Subsonic outflow
                    s = T / std::pow(rho, gamInf - 1.0);
                } else {  // Subsonic inflow
                    s = TInf / std::pow(rhoInf, gamInf - 1.0);
                    ut = utInf;
                    vt = vtInf;
                }

                un = 0.5 * (R1 + R2);
                c = 0.25 * (gamInf - 1.0) * (R1 - R2);
                T = c * c / gamInf;
                rho = std::pow(T / s, 1.0 / (gamInf - 1.0));

                // Recompute velocity components
                u = un * nx + ut;
                v = un * ny + vt;

                qB[0](j, k) = rho;
                qB[1](j, k) = rho * u;
                qB[2](j, k) = rho * v;
                qB[3](j, k) = rho * (T / (gamInf - 1.0) + 0.5 * (u * u + v * v));
            }
        }
    }

    return qB;
}
