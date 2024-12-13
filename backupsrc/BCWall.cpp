#include <Eigen/Dense>
#include <vector>
#include <cmath>

// Use Eigen for matrices and vectors
using namespace Eigen;
typedef MatrixXd Matrix;
typedef std::vector<Matrix> Matrix3D; // 3D vector for qB (4 matrices for 4 state variables)

// Function: bcwall
Matrix3D bcwall(const Matrix3D &q, const VectorXd &qInfPrim, const Matrix &xx, const Matrix &xy, int NJ, int NK) {
    // Initialize output boundary condition matrix qB
    Matrix3D qB(4, Matrix::Zero(NJ, NK));

    // Free-stream quantities
    double rhoInf = qInfPrim(0);
    double uInf = qInfPrim(1);
    double vInf = qInfPrim(2);
    double TInf = qInfPrim(3);
    double gamInf = qInfPrim(4);

    double h0Inf = TInf * gamInf / (gamInf - 1.0) + 0.5 * (uInf * uInf + vInf * vInf);

    // Loop over grid points
    for (int j = 0; j < NJ; ++j) {
        for (int k = 0; k < NK; ++k) {
            double rho = q[0](j, k);
            double u = q[1](j, k) / rho;
            double v = q[2](j, k) / rho;
            double T = (gamInf - 1.0) * (q[3](j, k) / rho - 0.5 * (u * u + v * v));
            double h0 = q[3](j, k) / rho + T;
            double s = T / std::pow(rho, gamInf - 1.0);
            double P = rho * T;

            // Compute normal and tangential components of velocity
            double normFactor = std::sqrt(xx(j, k) * xx(j, k) + xy(j, k) * xy(j, k));
            double nx = xx(j, k) / normFactor;
            double ny = xy(j, k) / normFactor;

            double un = u * nx + v * ny;  // Normal velocity
            double ut = u - un * nx;      // Tangential velocity in x
            double vt = v - un * ny;      // Tangential velocity in y

            // Set temperature based on total enthalpy and entropy
            T = (h0 - 0.5 * (ut * ut + vt * vt)) * (gamInf - 1.0) / gamInf;
            rho = std::pow(T / s, 1.0 / (gamInf - 1.0));

            // Update boundary condition state
            qB[0](j, k) = rho;
            qB[1](j, k) = rho * ut;  // Tangential momentum in x
            qB[2](j, k) = rho * vt;  // Tangential momentum in y
            qB[3](j, k) = rho * (T / (gamInf - 1.0) + 0.5 * (ut * ut + vt * vt));
        }
    }

    return qB;
}
