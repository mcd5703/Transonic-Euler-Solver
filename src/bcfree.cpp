#include "boundaryconditions.hpp"
#include <cmath>
#include <vector>

std::vector<std::vector<std::vector<double>>> bcfree(
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

    // cInf: Since TInf = 1/gamInf, cInf = sqrt(gamInf*TInf) = sqrt(gamInf*(1/gamInf))=1.
    // Matches MATLAB assumption cInf=1.

    std::vector<std::vector<std::vector<double>>> qB(NJ,
        std::vector<std::vector<double>>(NK,std::vector<double>(4,0.0)));

    for (int j = 0; j < NJ; j++) {
        for (int k = 0; k < NK; k++) {
            double rho = q[j][k][0];
            double rhou= q[j][k][1];
            double rhov= q[j][k][2];
            double rhoE= q[j][k][3];

            double u = rhou/rho;
            double v = rhov/rho;
            double T = (gamInf - 1.0)*((rhoE/rho) - 0.5*(u*u + v*v));

            double denom = std::sqrt(xx[j][k]*xx[j][k] + xy[j][k]*xy[j][k]);
            double nx = xx[j][k]/denom;
            double ny = xy[j][k]/denom;

            double un = u*nx + v*ny;
            double unInf = uInf*nx + vInf*ny;
            double c = std::sqrt(gamInf*T);
            double cInf = 1.0;

            double ut = u - un*nx;
            double vt = v - un*ny;
            double utInf = uInf - unInf*nx;
            double vtInf = vInf - unInf*ny;

            if (unInf + cInf < 0) {
                // Supersonic inflow
                qB[j][k][0] = rhoInf;
                qB[j][k][1] = rhoInf*uInf;
                qB[j][k][2] = rhoInf*vInf;
                qB[j][k][3] = rhoInf*(TInf/(gamInf-1.0) + 0.5*(uInf*uInf + vInf*vInf));
            } else if (un > c) {
                // Supersonic outflow - just copy q
                qB[j][k][0] = q[j][k][0];
                qB[j][k][1] = q[j][k][1];
                qB[j][k][2] = q[j][k][2];
                qB[j][k][3] = q[j][k][3];
            } else {
                // Subsonic
                double s = T/std::pow(rho,(gamInf-1.0));

                if (un > 0) {
                    // subsonic outflow: s is local, do nothing special
                    // s stays as is
                } else {
                    // subsonic inflow: use freestream s
                    s = TInf/std::pow(rhoInf,(gamInf-1.0));
                    ut = utInf;
                    vt = vtInf;
                }

                double R1 = un + 2*c/(gamInf-1.0);
                double R2 = unInf - 2*cInf/(gamInf-1.0);

                double unNew = 0.5*(R1 + R2);
                double cNew = 0.25*(gamInf - 1.0)*(R1 - R2);
                double Tnew = (cNew*cNew)/gamInf;
                double rhoNew = std::pow(Tnew/s, 1.0/(gamInf-1.0));
                double uNew = unNew*nx + ut;
                double vNew = unNew*ny + vt;

                qB[j][k][0] = rhoNew;
                qB[j][k][1] = rhoNew*uNew;
                qB[j][k][2] = rhoNew*vNew;
                qB[j][k][3] = rhoNew*(Tnew/(gamInf-1.0) + 0.5*(uNew*uNew + vNew*vNew));
            }
        }
    }

    return qB;
}
