#include "boundaryconditions.hpp"
#include <cmath>
#include <vector>

std::vector<std::vector<std::vector<double>>> bcwall(
    const std::vector<std::vector<std::vector<double>>> &q,
    const std::vector<double> &qInfPrim,
    const std::vector<std::vector<double>> &xx,
    const std::vector<std::vector<double>> &xy,
    const std::vector<std::vector<double>> &vol,
    int NJ, int NK)
{
    double rhoInf = qInfPrim[0];
    double uInf   = qInfPrim[1];
    double vInf   = qInfPrim[2];
    double TInf   = qInfPrim[3];
    double gamInf = qInfPrim[4];

    double h0Inf = TInf*gamInf/(gamInf - 1.0) + 0.5*(uInf*uInf + vInf*vInf);

    std::vector<std::vector<std::vector<double>>> qB(NJ,
        std::vector<std::vector<double>>(NK,std::vector<double>(4,0.0)));

    for (int j = 0; j < NJ; j++) {
        for (int k = 0; k < NK; k++) {
            double rho = q[j][k][0];
            double u   = q[j][k][1]/rho;
            double v   = q[j][k][2]/rho;
            double T   = (gamInf - 1.0)*((q[j][k][3]/rho) - 0.5*(u*u + v*v));
            double h0  = (q[j][k][3]/rho) + T;
            double s   = T/std::pow(rho,(gamInf-1.0));
            double P   = rho*T;
            double denom = std::sqrt(xx[j][k]*xx[j][k] + xy[j][k]*xy[j][k]);
            double nx = xx[j][k]/denom;
            double ny = xy[j][k]/denom;
            double un = u*nx + v*ny;
            double ut = u - un*nx;
            double vt = v - un*ny;

            // Set temperature based on total enthalpy and entropy
            T = (h0 - 0.5*(ut*ut + vt*vt))*(gamInf - 1.0)/gamInf;
            rho = std::pow(T/s, 1.0/(gamInf - 1.0));

            qB[j][k][0] = rho;
            qB[j][k][1] = rho*ut;
            qB[j][k][2] = rho*vt;
            qB[j][k][3] = rho*(T/(gamInf-1.0) + 0.5*(ut*ut + vt*vt));
        }
    }

    return qB;
}
