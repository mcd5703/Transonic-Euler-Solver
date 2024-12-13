#include "grid.hpp"
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>

// This implementation closely follows the MATLAB code provided, with appropriate 
// indexing and complex number handling in C++.

// MATLAB references: (1-based indexing in MATLAB, converted to 0-based here)
// In MATLAB:
// j: 1 to NJ, k: 1 to NK
// In C++:
// j: 0 to NJ-1, k: 0 to NK-1

void makegrid(int NACA, double R2, int NJ, int NK,
                std::vector<std::vector<double>> &x,
                std::vector<std::vector<double>> &y)
{
    // Extract NACA parameters
    int d1 = NACA / 1000;
    double m = d1 / 100.0;
    NACA = NACA - 1000 * d1;
    int d2 = NACA / 100;
    double p = d2 / 10.0;
    NACA = NACA - 100 * d2;
    double t = NACA / 100.0;

    double eps = 4 * std::sqrt(3.0) * t / 9.0;

    double z1 = 1.0;
    double denom = ((1 + 2 * eps) * (1 + 2 * eps)) + 2 * (1 + 2 * eps) + 1;
    double z2 = z1 - 4 * (1 + 2 * eps) / denom;

    // Set up target phi space
    std::vector<double> phi(NK,0.0);
    for (int k=0;k<NK;k++){
        phi[k] = 2*M_PI*( (double)k / (double)(NK-1) );
    }

    // Initial guess at theta
    std::vector<double> theta = phi;
    double errorMax = 999.0;

    // We'll store temporary arrays for the top row of x,y during iteration
    // Because MATLAB sets x(1,k),y(1,k), in C++ we use x[0][k], y[0][k].
    // We'll also need complex numbers for mapping
    typedef std::complex<double> dcomp;
    dcomp Z1(z1,0.0), Z2(z2,0.0);

    while (errorMax > 1e-10) {
        // Construct the airfoil surface on the top row j=0
        for (int k=0; k<NK; k++){
            x[0][k] = 0.5*(1 + std::cos(theta[k]));
            double val = x[0][k];
            double yt = 5*t*(0.2969*std::sqrt(val) - 0.1260*val - 0.3516*val*val 
                             + 0.2843*val*val*val - 0.1036*val*val*val*val);
            if (theta[k]>M_PI) yt = -yt;

            // Add camber line correction
            if (val < p) {
                y[0][k] = yt + m*(2*p*val - val*val)/(p*p);
            } else {
                y[0][k] = yt + m*((1 - 2*p) + 2*p*val - val*val)/((1 - p)*(1 - p));
            }
        }

        // Force endpoints to have y=0 as in MATLAB
        y[0][0] = 0.0;
        y[0][NK-1] = 0.0;

        // Trailing-edge angle
        double thetaU = std::atan2(y[0][1]-y[0][0], x[0][1]-x[0][0]);
        double thetaL = std::atan2(y[0][NK-2]-y[0][NK-1], x[0][NK-2]-x[0][NK-1]);
        if (thetaU < 0) thetaU += 2*M_PI;
        if (thetaL < 0) thetaL += 2*M_PI;
        double tauTE = thetaL - thetaU;
        double PTE = M_PI/(2*M_PI - tauTE);

        // Map the surface
        std::vector<std::vector<dcomp>> zeta(NJ, std::vector<dcomp>(NK,dcomp(0,0)));
        dcomp zeta1(-0.5*(z2 - z1)*PTE,0.0);
        dcomp zeta2(0.5*(z2 - z1)*PTE,0.0);

        for (int k=0; k<NK; k++){
            dcomp Z(x[0][k], y[0][k]);
            dcomp RHS = std::pow((Z - Z1)/(Z - Z2), PTE);

            zeta[0][k] = (zeta1 - RHS*zeta2)/(dcomp(1,0)-RHS);
            if ((theta[k]>M_PI) && (std::imag(zeta[0][k])>0)) {
                // Branch fix
                dcomp factor = std::exp(dcomp(0,-2.0*M_PI*PTE));
                RHS = RHS*factor;
                zeta[0][k] = (zeta1 - RHS*zeta2)/(dcomp(1,0)-RHS);
            }
        }

        // Find origin (zeta0)
        dcomp zeta0(0,0);
        for (int kk=0; kk<NK-1; kk++){
            double dth = theta[kk+1]-theta[kk];
            zeta0 += (zeta[0][kk+1]+zeta[0][kk])*dth;
        }
        zeta0 = zeta0/(2*M_PI);

        std::vector<double> phiA(NK,0.0);
        for (int kk=0; kk<NK; kk++){
            dcomp diff = zeta[0][kk]-zeta0;
            phiA[kk] = std::atan2(std::imag(diff), std::real(diff));
            if ((theta[kk]-M_PI/2)>0 && phiA[kk]<0) phiA[kk]+=2*M_PI;
        }

        double phi0 = phiA[0];
        phiA[NK-1] = phi0 + 2*M_PI;

        // error = phiA - (phi + phi0)
        std::vector<double> error(NK,0.0);
        for (int kk=0; kk<NK; kk++){
            error[kk] = phiA[kk] - (phi[kk] + phi0);
        }
        errorMax=0.0;
        for (int kk=0; kk<NK; kk++){
            double val = std::fabs(error[kk]);
            if(val>errorMax) errorMax=val;
        }

        // sensitivities and update
        std::vector<double> dthetaArr(NK,0.0);
        for (int kk=1; kk<NK-1; kk++){
            double dphidtheta = (phiA[kk+1]-phiA[kk-1])/(theta[kk+1]-theta[kk-1]);
            dthetaArr[kk] = -error[kk]/dphidtheta;
        }

        for (int kk=0; kk<NK; kk++){
            theta[kk]+= dthetaArr[kk];
        }
    }

    // Build pseudo-circular grid
    // We have final theta and phiA from the last iteration, we must reconstruct them:
    // Actually, at this point MATLAB code re-uses q computed in iteration. We can replicate final steps:
    // We'll re-do the surface mapping steps to get final zetaRow and phiA

    // Recompute surface with final theta
    for (int k=0; k<NK; k++){
        x[0][k] = 0.5*(1 + std::cos(theta[k]));
        double val = x[0][k];
        double yt = 5*t*(0.2969*std::sqrt(val)-0.1260*val -0.3516*val*val 
                         +0.2843*val*val*val -0.1036*val*val*val*val);
        if(theta[k]>M_PI) yt=-yt;
        double yc;
        if(val<p) {
            yc=m*(2*p*val - val*val)/(p*p);
        } else {
            yc=m*((1-2*p)+2*p*val - val*val)/((1-p)*(1-p));
        }
        y[0][k]=yt+yc;
        if(k==0) y[0][k]=0;
        if(k==NK-1) y[0][k]=0;
    }

    double thetaU = std::atan2(y[0][1]-y[0][0], x[0][1]-x[0][0]);
    double thetaL = std::atan2(y[0][NK-2]-y[0][NK-1], x[0][NK-2]-x[0][NK-1]);
    if(thetaU<0) thetaU+=2*M_PI;
    if(thetaL<0) thetaL+=2*M_PI;
    double tauTE = thetaL -thetaU;
    double PTE = M_PI/(2*M_PI - tauTE);

    dcomp zeta1(-0.5*(z2 - z1)*PTE,0.0);
    dcomp zeta2(0.5*(z2 - z1)*PTE,0.0);

    std::vector<dcomp> zetaRow(NK,dcomp(0,0));
    for (int k=0;k<NK;k++){
        dcomp Z(x[0][k],y[0][k]);
        dcomp RHS = std::pow((Z - Z1)/(Z - Z2),PTE);
        dcomp zeta_candidate = (zeta1 - RHS*zeta2)/(dcomp(1,0)-RHS);
        if ((theta[k]>M_PI)&&(std::imag(zeta_candidate)>0)) {
            dcomp factor = std::exp(dcomp(0,-2.0*M_PI*PTE));
            RHS = RHS*factor;
            zeta_candidate = (zeta1 - RHS*zeta2)/(dcomp(1,0)-RHS);
        }
        zetaRow[k]=zeta_candidate;
    }

    dcomp zeta0(0,0);
    for (int k=0;k<NK-1;k++){
        double dth = theta[k+1]-theta[k];
        zeta0 += (zetaRow[k+1]+zetaRow[k])*dth;
    }
    zeta0 = zeta0/(2*M_PI);

    std::vector<double> phiA(NK,0.0);
    for (int k=0;k<NK;k++){
        dcomp diff = zetaRow[k]-zeta0;
        phiA[k]=std::atan2(std::imag(diff),std::real(diff));
        if((theta[k]-M_PI/2)>0 && phiA[k]<0) phiA[k]+=2*M_PI;
    }

    // Now create the entire (NJ x NK) grid
    // j goes 0 to NJ-1, k goes 0 to NK-1
    // rad = R1*R2/(R2 + ((j)/(NJ-1))*(R1-R2))
    // zeta(j,k) = zeta0 + rad*cos(phiA(k)) + i*rad*sin(phiA(k))
    // Map back

    for (int j=0;j<NJ;j++){
        for (int k=0;k<NK;k++){
            double R1 = std::abs(zetaRow[k]-zeta0);
            double rad = R1*R2/(R2 + ((double)j/(NJ-1))*(R1 - R2));
            dcomp zetaP = zeta0 + rad*std::cos(phiA[k]) + dcomp(0,1)*rad*std::sin(phiA[k]);

            dcomp RHS = std::pow((zetaP - zeta1)/(zetaP - zeta2),1.0/PTE);
            dcomp ZZ = (Z2*RHS - Z1)/(RHS - dcomp(1,0));

            x[j][k]=std::real(ZZ);
            y[j][k]=std::imag(ZZ);
        }
    }
}
