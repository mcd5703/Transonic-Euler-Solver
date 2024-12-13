#ifndef RESID_HPP
#define RESID_HPP

#include <vector>
#include <cmath>
#include <algorithm>

// Include the header that provides HLLCFlux function
#include "HLLCFlux.hpp"

inline void resid(
    const std::vector<std::vector<std::vector<double>>> &q,
    const std::vector<std::vector<double>> &xx,
    const std::vector<std::vector<double>> &xy,
    const std::vector<std::vector<double>> &yx,
    const std::vector<std::vector<double>> &yy,
    const std::vector<std::vector<double>> &vol,
    const std::vector<double> &qInfPrim,
    int NJ, int NK, double CFL,
    std::vector<std::vector<std::vector<double>>> &R,
    std::vector<std::vector<double>> &dt
) {
    // Extract freestream properties
    double rhoInf = qInfPrim[0];
    double uInf = qInfPrim[1];
    double vInf = qInfPrim[2];
    double TInf = qInfPrim[3];
    double gamInf = qInfPrim[4];
    double h0Inf = TInf * gamInf / (gamInf - 1.) + 0.5 * (uInf * uInf + vInf * vInf);

    // Initialize R and dt
    for (int j = 0; j < NJ; j++) {
        for (int k = 0; k < NK; k++) {
            for (int n = 0; n < 4; n++) {
                R[j][k][n] = 0.0;
            }
            dt[j][k] = 0.0;
        }
    }

    std::vector<std::vector<double>> sigj(NJ, std::vector<double>(NK, 0.0));
    std::vector<std::vector<double>> sigk(NJ, std::vector<double>(NK, 0.0));

    // Compute spectral radii
    for (int j = 1; j < NJ - 1; j++) {
        for (int k = 0; k < NK; k++) {
            double rho = q[j][k][0];
            double u = q[j][k][1] / rho;
            double v = q[j][k][2] / rho;
            double T = (gamInf - 1.) * ((q[j][k][3] / rho) - 0.5 * (u * u + v * v));
            sigj[j][k] = std::fabs(xx[j][k]*u + xy[j][k]*v) +
                         std::sqrt((xx[j][k]*xx[j][k] + xy[j][k]*xy[j][k]) * 1.4 * T);
            sigk[j][k] = std::fabs(yx[j][k]*u + yy[j][k]*v) +
                         std::sqrt((yx[j][k]*yx[j][k] + yy[j][k]*yy[j][k]) * 1.4 * T);
            double CFLeff = std::max(CFL/(1 + std::sqrt(vol[j][k])), 1.0);
            dt[j][k] = CFLeff * vol[j][k] / (std::max(sigj[j][k], sigk[j][k]));
        }
    }

    // Helper lambda for limiter
    auto limiter = [&](double a, double b){
        if (b == 0.0) return 1e8;
        return (a/b);
    };

    // HLLC wrapper lambda
    auto HLLCfun = [&](const std::vector<double>& qa, const std::vector<double>& qb, double xxA, double xyA, double volA){
        std::vector<double> F(4,0.0);
        HLLCFlux(qa, qb, gamInf, xxA, xyA, volA, F);
        // subtract free-stream flux
        double UU = xxA*uInf + xyA*vInf;
        F[0] -= rhoInf*UU;
        F[1] -= rhoInf*UU*uInf + xxA*rhoInf*TInf;
        F[2] -= rhoInf*UU*vInf + xyA*rhoInf*TInf;
        F[3] -= rhoInf*UU*h0Inf;
        return F;
    };

    // Compute K residuals (Periodic in K)
    for (int j = 1; j < NJ - 1; j++) {
        for (int k = 0; k < NK - 1; k++) {
            double yxA = 0.5 * (yx[j][k] + yx[j][k+1]);
            double yyA = 0.5 * (yy[j][k] + yy[j][k+1]);
            double volA = 0.5 * (vol[j][k] + vol[j][k+1]);

            int km = (k == 0) ? (NK - 2) : (k - 1);
            int kp = k + 1;
            int kpp = (k + 2 >= NK) ? 1 : k + 2;

            std::vector<double> qL(4), qR(4);
            for (int n = 0; n < 4; n++) {
                double r = limiter(q[j][kp][n] - q[j][k][n], q[j][k][n] - q[j][km][n]);
                double phi = std::max(0.0, std::min({2*r, (1+2*r)/3.0, 2.0}));
                qL[n] = q[j][k][n] + 0.5 * phi * (q[j][k][n]-q[j][km][n]);

                r = limiter(q[j][k][n] - q[j][kp][n], q[j][kp][n]-q[j][kpp][n]);
                phi = std::max(0.0, std::min({2*r, (1+2*r)/3.0, 2.0}));
                qR[n] = q[j][kp][n] + 0.5 * phi * (q[j][kp][n]-q[j][kpp][n]);
            }

            auto fhat = HLLCfun(qL,qR,yxA,yyA,volA);
            for (int n = 0; n < 4; n++) {
                R[j][k][n]   += fhat[n];
                R[j][k+1][n] -= fhat[n];
            }
        }
        for (int n = 0; n < 4; n++) {
            R[j][0][n] += R[j][NK-1][n];
            R[j][NK-1][n] = R[j][0][n];
        }
    }

    // Compute J residuals
    for (int k = 0; k < NK; k++) {
        for (int j = 0; j < NJ - 1; j++) {
            double xxA = 0.5 * (xx[j][k] + xx[j+1][k]);
            double xyA = 0.5 * (xy[j][k] + xy[j+1][k]);
            double volA = 0.5 * (vol[j][k] + vol[j+1][k]);

            int jm = (j == 0) ? 0 : j - 1;
            int jp = j + 1;
            int jpp = (j + 2 >= NJ) ? NJ - 1 : j + 2;

            std::vector<double> qL(4), qR(4);
            for (int n = 0; n < 4; n++) {
                double r = limiter(q[jp][k][n] - q[j][k][n], q[j][k][n] - q[jm][k][n]);
                double phi = std::max(0.0, std::min({2*r, (1+2*r)/3.0, 2.0}));
                qL[n] = q[j][k][n] + 0.5 * phi * (q[j][k][n]-q[jm][k][n]);

                r = limiter(q[j][k][n] - q[jp][k][n], q[jp][k][n]-q[jpp][k][n]);
                phi = std::max(0.0, std::min({2*r, (1+2*r)/3.0, 2.0}));
                qR[n] = q[jp][k][n] + 0.5 * phi * (q[jp][k][n]-q[jpp][k][n]);
            }

            auto fhat = HLLCfun(qL,qR,xxA,xyA,volA);
            for (int n = 0; n < 4; n++) {
                R[j][k][n]   += fhat[n];
                R[j+1][k][n] -= fhat[n];
            }
        }
    }

    // Volume scale and multiply by dt
    for (int k = 0; k < NK; k++) {
        for (int j = 1; j < NJ - 1; j++) {
            for (int n = 0; n < 4; n++) {
                R[j][k][n] = dt[j][k] * R[j][k][n] / vol[j][k];
            }
        }
        for (int n = 0; n < 4; n++) {
            R[0][k][n] = 0;
            R[NJ-1][k][n] = 0;
        }
    }
}

#endif // RESID_HPP
