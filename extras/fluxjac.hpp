#ifndef FLUXJAC_HPP
#define FLUXJAC_HPP

#include <vector>
#include <cmath>

inline void fluxjac(
    const std::vector<double> &v1,
    const std::vector<double> &q,
    double xx, double xy, double gamInf,
    std::vector<double> &fout
) {
    // q = [rho, rho*u, rho*v, rho*E]
    double rho = q[0];
    double u = q[1]/rho;
    double v = q[2]/rho;
    double e0 = q[3]/rho;  // e0 = total energy per unit mass
    double phi2 = 0.5*(gamInf - 1)*(u*u + v*v);
    double UU = xx*u + xy*v;

    // Construct matrix A based on the given formula
    // A(1,:) = [0, xx, xy, 0]
    // A(2,:) = [xx*phi2 - UU*u, UU - xx*(gamInf-2)*u, xy*u - (gamInf-1)*xx*v, xx*(gamInf-1)]
    // A(3,:) = [xy*phi2 - UU*v, xx*v - xy*(gamInf-1)*u, UU - xy*(gamInf-2)*v, xy*(gamInf-1)]
    // A(4,:) = [UU*(2*phi2 - gamInf*e0), xx*(gamInf*e0 - phi2) - (gamInf-1)*u*UU, xy*(gamInf*e0 - phi2) - (gamInf-1)*v*UU, gamInf*UU]

    double A[4][4];
    A[0][0] = 0.0; 
    A[0][1] = xx; 
    A[0][2] = xy; 
    A[0][3] = 0.0;

    A[1][0] = xx*phi2 - UU*u;
    A[1][1] = UU - xx*(gamInf-2)*u;
    A[1][2] = xy*u - (gamInf-1)*xx*v;
    A[1][3] = xx*(gamInf-1);

    A[2][0] = xy*phi2 - UU*v;
    A[2][1] = xx*v - xy*(gamInf-1)*u;
    A[2][2] = UU - xy*(gamInf-2)*v;
    A[2][3] = xy*(gamInf-1);

    A[3][0] = UU*(2*phi2 - gamInf*e0);
    A[3][1] = xx*(gamInf*e0 - phi2) - (gamInf-1)*u*UU;
    A[3][2] = xy*(gamInf*e0 - phi2) - (gamInf-1)*v*UU;
    A[3][3] = gamInf*UU;

    // Compute A*v1
    fout.assign(4,0.0);
    for (int n=0; n<4; n++) {
        for (int k=0; k<4; k++) {
            fout[n] += A[n][k]*v1[k];
        }
    }
}

#endif // FLUXJAC_HPP
