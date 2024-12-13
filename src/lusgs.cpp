#include "lusgs.hpp"
#include <cmath>
#include <vector>
#include <algorithm>
#include "simpleFunctions.hpp" // for fluxjac
// qdotq not needed here directly, but fluxjac is needed.

void lusgs(
    const std::vector<std::vector<std::vector<double>>> &R,
    const std::vector<std::vector<std::vector<double>>> &q,
    const std::vector<std::vector<double>> &xx,
    const std::vector<std::vector<double>> &xy,
    const std::vector<std::vector<double>> &yx,
    const std::vector<std::vector<double>> &yy,
    const std::vector<std::vector<double>> &vol,
    const std::vector<std::vector<double>> &dt,
    double gamInf, int NJ, int NK,
    std::vector<std::vector<std::vector<double>>> &dq
)
{
    // dq and dqs are [NJ][NK][4]
    // Initialize dq and dqs to zero
    std::vector<std::vector<std::vector<double>>> dqs(NJ, std::vector<std::vector<double>>(NK,std::vector<double>(4,0.0)));

    std::vector<std::vector<double>> sigj(NJ, std::vector<double>(NK,0.0));
    std::vector<std::vector<double>> sigk(NJ, std::vector<double>(NK,0.0));

    // Build spectral radii
    for (int j=0; j<NJ; j++) {
        for (int k=0; k<NK; k++) {
            double rho = q[j][k][0];
            double u = q[j][k][1]/rho;
            double v = q[j][k][2]/rho;
            double T = (gamInf - 1.0)*( (q[j][k][3]/rho) - 0.5*(u*u + v*v));
            sigj[j][k] = std::fabs(xx[j][k]*u + xy[j][k]*v) + std::sqrt((xx[j][k]*xx[j][k] + xy[j][k]*xy[j][k])*1.4*T);
            sigk[j][k] = std::fabs(yx[j][k]*u + yy[j][k]*v) + std::sqrt((yx[j][k]*yx[j][k] + yy[j][k]*yy[j][k])*1.4*T);
            // Now divide by vol(j,k)
            sigj[j][k] = sigj[j][k]/vol[j][k];
            sigk[j][k] = sigk[j][k]/vol[j][k];
        }
    }

    // Temporary vectors for fluxjac
    std::vector<double> tmp1(4,0.0);
    std::vector<double> tmp2(4,0.0);

    // Forward substitution loop
    for (int j=0; j<NJ; j++) {
        for (int k=0; k<NK; k++) {
            int jm = j-1;
            int km = k-1;
            double facx = 2.0;
            double facy = 2.0;

            // if (j > 1) in MATLAB means j>0 in 0-based
            if (j > 0) {
                facx = 1.0;
                // fluxjac(dqs(jm,k,:),q(jm,k,:),xx(jm,k)/vol(jm,k),xy(jm,k)/vol(jm,k),gamInf)
                // Extract dqs(jm,k,:) and q(jm,k,:)
                std::vector<double> v1(4), qq(4);
                for (int n=0;n<4;n++){
                    v1[n]=dqs[jm][k][n];
                    qq[n]=q[jm][k][n];
                }
                fluxjac(v1, qq, xx[jm][k]/vol[jm][k], xy[jm][k]/vol[jm][k], gamInf, tmp1);
                for (int n=0;n<4;n++){
                    tmp1[n] = 0.5*(tmp1[n] + sigj[jm][k]*dqs[jm][k][n]);
                }
            } else {
                for (int n=0;n<4;n++) tmp1[n]=0.0;
            }

            if (k > 0) {
                facy = 1.0;
                std::vector<double> v1(4), qq(4);
                for (int n=0;n<4;n++){
                    v1[n]=dqs[j][km][n];
                    qq[n]=q[j][km][n];
                }
                fluxjac(v1, qq, yx[j][km]/vol[j][km], yy[j][km]/vol[j][km], gamInf, tmp2);
                for (int n=0;n<4;n++){
                    tmp2[n] = 0.5*(tmp2[n] + sigk[j][km]*dqs[j][km][n]);
                }
            } else {
                for (int n=0;n<4;n++) tmp2[n]=0.0;
            }

            double D = 1.0 + dt[j][k]*(facx*sigj[j][k] + facy*sigk[j][k]);

            for (int n=0;n<4;n++){
                dqs[j][k][n] = (-R[j][k][n]*D*vol[j][k] + dt[j][k]*(tmp1[n]+tmp2[n]))/D;
            }
        }
    }

    // Clean up periodicity in dqs
    for (int j=0;j<NJ;j++){
        for (int n=0;n<4;n++){
            dqs[j][0][n] = 0.5*(dqs[j][0][n] + dqs[j][NK-1][n]);
            dqs[j][NK-1][n] = dqs[j][0][n];
        }
    }

    // Backward substitution loop
    for (int j=NJ-1; j>=0; j--) {
        for (int k=NK-1; k>=0; k--) {
            int jp = j+1;
            int kp = k+1;
            double facx = 2.0;
            double facy = 2.0;

            // if (j < NJ) in MATLAB means j<NJ => j<=NJ-1 always true for j in range.
            // Actually it meant if not the last row. In MATLAB indexing j runs from 1:Nj.
            // The condition allows flux calculation with jp. 
            // If j<NJ in MATLAB means if j <= NJ-1 in 1-based indexing means j<NJ (no eq).
            // In C++ j runs from NJ-1 down, 
            // j<NJ is always true. We must interpret original logic:
            // The code tries to compute tmp1 using dq(jp,k,...) if j<NJ means j<=(NJ-1) since max j=NJ.
            // This is always true except when j=NJ in MATLAB, which doesn't happen since loop stops at j=1.
            // So condition is always true except at last row? 
            // In C++ last row is j=NJ-1. If jp=j+1 must be valid: if j=NJ-1, jp=NJ out of range.
            // So condition (j < NJ) in MATLAB means (j < NJ) 1-based => j<=NJ-1. For safety we use (j<NJ-1)
            if (j < NJ-1) {
                facx=1.0;
                std::vector<double> v1(4), qq(4);
                for (int n=0;n<4;n++){
                    v1[n]=dq[jp][k][n];
                    qq[n]=q[jp][k][n];
                }
                fluxjac(v1, qq, xx[jp][k]/vol[jp][k], xy[jp][k]/vol[jp][k], gamInf, tmp1);
                for (int n=0;n<4;n++){
                    tmp1[n]=0.5*(tmp1[n] - sigj[jp][k]*dq[jp][k][n]);
                }
            } else {
                for (int n=0;n<4;n++) tmp1[n]=0.0;
            }

            if (k < NK-1) {
                facy=1.0;
                std::vector<double> v1(4), qq(4);
                for (int n=0;n<4;n++){
                    v1[n]=dq[j][kp][n];
                    qq[n]=q[j][kp][n];
                }
                fluxjac(v1, qq, yx[j][kp]/vol[j][kp], yy[j][kp]/vol[j][kp], gamInf, tmp2);
                for (int n=0;n<4;n++){
                    tmp2[n]=0.5*(tmp2[n] - sigk[j][kp]*dq[j][kp][n]);
                }
            } else {
                for (int n=0;n<4;n++) tmp2[n]=0.0;
            }

            double D = 1.0 + dt[j][k]*(facx*sigj[j][k] + facy*sigk[j][k]);

            for (int n=0;n<4;n++){
                dq[j][k][n] = (dqs[j][k][n] - dt[j][k]*(tmp1[n]+tmp2[n]))/D;
            }
        }
    }

    // Clean up periodicity in dq
    for (int j=0;j<NJ;j++){
        for (int n=0;n<4;n++){
            dq[j][0][n] = 0.5*(dq[j][0][n] + dq[j][NK-1][n]);
            dq[j][NK-1][n] = dq[j][0][n];
        }
    }

    // Divide dq by vol
    for (int j=0;j<NJ;j++){
        for (int k=0;k<NK;k++){
            for (int n=0;n<4;n++){
                dq[j][k][n] = dq[j][k][n]/vol[j][k];
            }
        }
    }
}
