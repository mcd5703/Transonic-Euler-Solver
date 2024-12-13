#include "resid.hpp"
#include <cmath>
#include <algorithm>
#include <vector>
#include "HLLCFlux.hpp" 

void resid(
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
)
{
    // Extract freestream properties
    double rhoInf = qInfPrim[0];
    double uInf = qInfPrim[1];
    double vInf = qInfPrim[2];
    double TInf = qInfPrim[3];
    double gamInf = qInfPrim[4];
    double h0Inf = TInf*gamInf/(gamInf-1.0) + 0.5*(uInf*uInf + vInf*vInf);

    // Initialize R and dt
    for (int j=0; j<NJ; j++){
        for (int k=0; k<NK; k++){
            for (int n=0; n<4; n++){
                R[j][k][n]=0.0;
            }
            dt[j][k]=0.0;
        }
    }

    std::vector<std::vector<double>> sigj(NJ,std::vector<double>(NK,0.0));
    std::vector<std::vector<double>> sigk(NJ,std::vector<double>(NK,0.0));

    // Build spectral radii
    for (int j=1; j<NJ-1; j++){
        for (int k=0; k<NK; k++){
            double rho = q[j][k][0];
            double u = q[j][k][1]/rho;
            double v = q[j][k][2]/rho;
            double T = (gamInf - 1.0)*((q[j][k][3]/rho)-0.5*(u*u+v*v));
            sigj[j][k] = std::fabs(xx[j][k]*u + xy[j][k]*v) + std::sqrt((xx[j][k]*xx[j][k] + xy[j][k]*xy[j][k])*1.4*T);
            sigk[j][k] = std::fabs(yx[j][k]*u + yy[j][k]*v) + std::sqrt((yx[j][k]*yx[j][k] + yy[j][k]*yy[j][k])*1.4*T);
            double CFLeff = std::max(CFL/(1+std::sqrt(vol[j][k])),1.0);
            dt[j][k] = CFLeff*vol[j][k]/(std::max(sigj[j][k],sigk[j][k]));
        }
    }

    // Use std::vector<double> for qL, qR, fhat to match HLLCFlux signature
    std::vector<double> qL(4), qR(4), fhat(4);

    // Do the K residuals (periodic)
    for (int j=1; j<NJ-1; j++){
        for (int k=0; k<NK-1; k++){
            // Periodic indexing in K direction
            int km = k-1;
            if (k == 0) km = NK-2; 
            int kp = k+1;
            int kpp = k+2;
            if (kpp >= NK) kpp = 1;

            double yxA = 0.5*(yx[j][k] + yx[j][k+1]);
            double yyA = 0.5*(yy[j][k] + yy[j][k+1]);
            double volA= 0.5*(vol[j][k] + vol[j][k+1]);

            // Left state
            for (int n=0; n<4; n++){
                double denomL = (q[j][k][n]-q[j][km][n]);
                double r = denomL==0 ? 1e8 : (q[j][kp][n] - q[j][k][n])/denomL;
                double phi = std::max(0.0,std::min({2*r,(1+2*r)/3.0,2.0}));
                qL[n] = q[j][k][n] + 0.5*phi*(q[j][k][n]-q[j][km][n]);
            }

            // Right state
            for (int n=0; n<4; n++){
                double denomR = (q[j][kp][n]-q[j][kpp][n]);
                double r = denomR==0 ? 1e8 : (q[j][k][n] - q[j][kp][n])/denomR;
                double phi = std::max(0.0,std::min({2*r,(1+2*r)/3.0,2.0}));
                qR[n] = q[j][kp][n] + 0.5*phi*(q[j][kp][n]-q[j][kpp][n]);
            }

            HLLCFlux(qL, qR, gamInf, yxA, yyA, volA, fhat);

            // Subtract free-stream fluxes
            double UU = yxA*uInf + yyA*vInf;
            fhat[0] -= rhoInf*UU;
            fhat[1] -= rhoInf*UU*uInf + yxA*rhoInf*TInf;
            fhat[2] -= rhoInf*UU*vInf + yyA*rhoInf*TInf;
            fhat[3] -= rhoInf*UU*h0Inf;

            for (int n=0;n<4;n++){
                R[j][k][n]   += fhat[n];
                R[j][k+1][n] -= fhat[n];
            }
        }
        for (int n=0;n<4;n++){
            R[j][0][n] += R[j][NK-1][n];
            R[j][NK-1][n] = R[j][0][n];
        }
    }

    // J residuals
    for (int k=0;k<NK;k++){
        for (int j=0; j<NJ-1; j++){
            int jm = j-1;
            if (j == 0) jm=0; 
            int jp = j+1;
            int jpp = j+2;
            if (jpp >= NJ) jpp=NJ-1;

            double xxA = 0.5*(xx[j][k] + xx[j+1][k]);
            double xyA = 0.5*(xy[j][k] + xy[j+1][k]);
            double volA= 0.5*(vol[j][k] + vol[j+1][k]);

            // Left state
            for (int n=0;n<4;n++){
                double denomL = (q[j][k][n]-q[jm][k][n]);
                double r = denomL==0 ? 1e8 : (q[jp][k][n]-q[j][k][n])/denomL;
                double phi = std::max(0.0,std::min({2*r,(1+2*r)/3.0,2.0}));
                qL[n] = q[j][k][n] + 0.5*phi*(q[j][k][n]-q[jm][k][n]);
            }

            // Right state
            for (int n=0;n<4;n++){
                double denomR = (q[jp][k][n]-q[jpp][k][n]);
                double r = denomR==0 ? 1e8 : (q[j][k][n]-q[jp][k][n])/denomR;
                double phi = std::max(0.0,std::min({2*r,(1+2*r)/3.0,2.0}));
                qR[n] = q[jp][k][n] + 0.5*phi*(q[jp][k][n]-q[jpp][k][n]);
            }

            HLLCFlux(qL, qR, gamInf, xxA, xyA, volA, fhat);

            // Subtract free-stream fluxes
            double UU = xxA*uInf + xyA*vInf;
            fhat[0] -= rhoInf*UU;
            fhat[1] -= rhoInf*UU*uInf + xxA*rhoInf*TInf;
            fhat[2] -= rhoInf*UU*vInf + xyA*rhoInf*TInf;
            fhat[3] -= rhoInf*UU*h0Inf;

            for (int n=0;n<4;n++){
                R[j][k][n]   += fhat[n];
                R[j+1][k][n] -= fhat[n];
            }
        }
    }

    // Volume-scale and multiply by DT
    for (int k=0; k<NK; k++){
        for (int j=1; j<NJ-1; j++){
            for (int n=0;n<4;n++){
                R[j][k][n] = dt[j][k]*R[j][k][n]/vol[j][k];
            }
        }
        for (int n=0;n<4;n++){
            R[0][k][n] = 0.0;
            R[NJ-1][k][n] = 0.0;
        }
    }
}
