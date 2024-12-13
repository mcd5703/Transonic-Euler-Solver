#include <iostream>
#include <cmath>
#include <vector>
#include <matplot/matplot.h> 
// Custom headers for our solver code
#include "solver.hpp"
// #include "grid.hpp"

void Solver(int NACA, double Minf, double alpha, double gamInf,
            int itmax, double CFL, int nDisp, int NJ, int NK, double R2)
{

    // Output the values for each variable in "config.h"
	std::cout << "-----------------------------------------------------------" << std::endl;
	std::cout << "NOTICE: The following defaults were read through main... " << std::endl;
    std::cout << "DEFAULT_NACA: " << NACA << std::endl;
    std::cout << "DEFAULT_MINF: " << Minf << std::endl;
    std::cout << "DEFAULT_ALPHA: " << alpha << std::endl;
    std::cout << "DEFAULT_GAMINF: " << gamInf << std::endl;
    std::cout << "DEFAULT_ITMAX: " << itmax << std::endl;
    std::cout << "DEFAULT_CFL: " << CFL << std::endl;
    std::cout << "DEFAULT_NDISP: " << nDisp << std::endl;
    std::cout << "DEFAULT_GRID_NJ: " << NJ << std::endl;
    std::cout << "DEFAULT_GRID_NK: " << NK << std::endl;
    std::cout << "DEFAULT_GRID_R2: " << R2 << std::endl;
	std::cout << "-----------------------------------------------------------" << std::endl;


    // // Make grid
    // std::vector<std::vector<double>> x(NJ, std::vector<double>(NK,0.0));
    // std::vector<std::vector<double>> y(NJ, std::vector<double>(NK,0.0));
    // makegrid(NACA,R2,NJ,NK,x,y);

    // // Compute metrics
    // std::vector<std::vector<double>> xx(NJ, std::vector<double>(NK,0.0));
    // std::vector<std::vector<double>> xy(NJ, std::vector<double>(NK,0.0));
    // std::vector<std::vector<double>> yx(NJ, std::vector<double>(NK,0.0));
    // std::vector<std::vector<double>> yy(NJ, std::vector<double>(NK,0.0));
    // std::vector<std::vector<double>> vol(NJ, std::vector<double>(NK,0.0));
    // metrics(x,y,NJ,NK,xx,xy,yx,yy,vol);

    // // Initialize solution to freestream
    // double rhoInf = 1.0;
    // double uInf = Minf*std::cos(alpha*M_PI/180.0);
    // double vInf = Minf*std::sin(alpha*M_PI/180.0);
    // double TInf = 1.0/gamInf;
    // std::vector<double> qInfPrim = {rhoInf, uInf, vInf, TInf, gamInf};

    // std::vector<std::vector<std::vector<double>>> q(NJ, std::vector<std::vector<double>>(NK,std::vector<double>(4,0.0)));
    // for (int j = 0; j < NJ; j++){
    //     for (int k = 0; k < NK; k++){
    //         q[j][k][0] = rhoInf; // Density
    //         q[j][k][1] = rhoInf*uInf; // XMomentum
    //         q[j][k][2] = rhoInf*vInf; // YMomentum
    //         q[j][k][3] = rhoInf*(TInf/(gamInf-1.0) + 0.5*(uInf*uInf + vInf*vInf));
    //     }
    // }

    // std::vector<std::vector<std::vector<double>>> R(NJ, std::vector<std::vector<double>>(NK,std::vector<double>(4,0.0)));
    // std::vector<std::vector<double>> dt(NJ, std::vector<double>(NK,0.0));
    // std::vector<std::vector<std::vector<double>>> qn = q;
    // std::vector<std::vector<std::vector<double>>> RN(NJ, std::vector<std::vector<double>>(NK,std::vector<double>(4,0.0)));
    // std::vector<std::vector<std::vector<double>>> v1(NJ, std::vector<std::vector<double>>(NK,std::vector<double>(4,0.0)));
    // std::vector<std::vector<std::vector<double>>> v2(NJ, std::vector<std::vector<double>>(NK,std::vector<double>(4,0.0)));

    // for (int it = 1; it <= itmax; it++){
    //     // q(1,:,:) = 1.5*q(2,:,:) - 0.5*q(3,:,:)
    //     for (int k = 0; k < NK; k++){
    //         for (int c = 0; c < 4; c++){
    //             q[0][k][c] = 1.5*q[1][k][c] - 0.5*q[2][k][c];
    //         }
    //     }

    //     bcwall(q, qInfPrim, xx[0], xy[0], vol[0], NK);

    //     // q(NJ,:,:) = q(NJ-1,:,:)
    //     for (int k = 0; k < NK; k++){
    //         for (int c = 0; c < 4; c++){
    //             q[NJ-1][k][c] = q[NJ-2][k][c];
    //         }
    //     }

    //     bcfree(q, qInfPrim, xx[NJ-1], xy[NJ-1], vol[NJ-1], NK);

    //     qn = q;

    //     resid(q,xx,xy,yx,yy,vol,qInfPrim,NJ,NK,CFL,R,dt);
    //     lusgs(R,q,xx,xy,yx,yy,vol,dt,gamInf,NJ,NK,v1);
    //     RN = R;

    //     double h11 = qdotq(v1,v1,NJ,NK);
    //     h11 = std::sqrt(h11);
    //     // v1 = v1/h11
    //     for (int j = 0; j < NJ; j++){
    //         for (int k = 0; k < NK; k++){
    //             for (int c = 0; c < 4; c++){
    //                 v1[j][k][c] /= h11;
    //             }
    //         }
    //     }

    //     double eps = 2e-8;
    //     // q = qn + (2e-8)*v1
    //     for (int j = 0; j < NJ; j++){
    //         for (int k = 0; k < NK; k++){
    //             for (int c = 0; c < 4; c++){
    //                 q[j][k][c] = qn[j][k][c] + eps*v1[j][k][c];
    //             }
    //         }
    //     }

    //     // Reapply BC
    //     for (int k = 0; k < NK; k++){
    //         for (int c = 0; c < 4; c++){
    //             q[0][k][c] = 1.5*q[1][k][c] - 0.5*q[2][k][c];
    //         }
    //     }

    //     bcwall(q, qInfPrim, xx[0], xy[0], vol[0], NK);

    //     for (int k = 0; k < NK; k++){
    //         for (int c = 0; c < 4; c++){
    //             q[NJ-1][k][c] = q[NJ-2][k][c];
    //         }
    //     }
    //     bcfree(q, qInfPrim, xx[NJ-1], xy[NJ-1], vol[NJ-1], NK);

    //     resid(q,xx,xy,yx,yy,vol,qInfPrim,NJ,NK,CFL,R,dt);

    //     // v2 = v1 + (R - RN)/(2e-8)
    //     for (int j = 0; j < NJ; j++){
    //         for (int k = 0; k < NK; k++){
    //             for (int c = 0; c < 4; c++){
    //                 v2[j][k][c] = v1[j][k][c] + (R[j][k][c] - RN[j][k][c])/eps;
    //             }
    //         }
    //     }

    //     double h22 = qdotq(v2,v2,NJ,NK);
    //     double h12 = qdotq(R,v2,NJ,NK);
    //     double a1 = -h12/h22;
    //     a1 = std::max(a1,0.5*h11);

    //     // q = qn + a1*v1
    //     for (int j = 0; j < NJ; j++){
    //         for (int k = 0; k < NK; k++){
    //             for (int c = 0; c < 4; c++){
    //                 q[j][k][c] = qn[j][k][c] + a1*v1[j][k][c];
    //             }
    //         }
    //     }

    //     if ((it % nDisp) == 0){
    //         // Compute residual norm
    //         double sumsq = 0.0;
    //         int count = 0;
    //         for (int j = 1; j < NJ-1; j++){
    //             for (int k = 0; k < NK; k++){
    //                 sumsq += R[j][k][0]*R[j][k][0];
    //                 count++;
    //             }
    //         }
    //         double rhoresid = std::sqrt(sumsq)/CFL;

    //         // Compute Mplot, P, Cp
    //         std::vector<std::vector<double>> Mplot(NJ,std::vector<double>(NK,0.0));
    //         std::vector<std::vector<double>> Pmat(NJ,std::vector<double>(NK,0.0));
    //         std::vector<std::vector<double>> Cp(NJ,std::vector<double>(NK,0.0));

    //         for (int j = 0; j < NJ; j++){
    //             for (int k = 0; k < NK; k++){
    //                 double rho = q[j][k][0];
    //                 double rhou = q[j][k][1];
    //                 double rhov = q[j][k][2];
    //                 double rhoE = q[j][k][3];
    //                 double u = rhou/rho;
    //                 double v = rhov/rho;
    //                 double p = (gamInf-1.0)*(rhoE - 0.5*(rhou*rhou + rhov*rhov)/rho);
    //                 double a2 = gamInf*(gamInf-1.0)*(rhoE/rho - 0.5*(u*u + v*v));
    //                 double M = std::sqrt((u*u+v*v)/a2);
    //                 Mplot[j][k] = M;
    //                 Pmat[j][k] = p;
    //                 Cp[j][k] = (p - 1.0/gamInf)/(0.5*Minf*Minf);
    //             }
    //         }

    //         // Integrate forces on j=1 -> j=0 in C++
    //         double cx = 0.0;
    //         double cy = 0.0;
    //         int J = 0;
    //         for (int K = 0; K < NK-1; K++){
    //             cx -= 0.5*(Cp[J][K]*xx[J][K] + Cp[J][K+1]*xx[J][K+1]);
    //             cy -= 0.5*(Cp[J][K]*xy[J][K] + Cp[J][K+1]*xy[J][K+1]);
    //         }

    //         double cl = cy*std::cos(alpha*M_PI/180.0) - cx*std::sin(alpha*M_PI/180.0);
    //         double cd = cx*std::cos(alpha*M_PI/180.0) + cy*std::sin(alpha*M_PI/180.0);

    //         std::cout << "Iter " << it << ": resid = " << rhoresid << ", cl = " << cl << ", cd = " << cd << std::endl;

    //         // Plotting using matplot++ (commented out)
    //         /*
    //         using namespace matplot;
    //         figure();
    //         // Convert Mplot and (x,y) into flat vectors for contourf if needed
    //         // ...
    //         // contourf(...) 
    //         // xlim({-1,2});
    //         // ylim({-1,1});
    //         // daspect({1,1,1});
    //         // drawnow();
    //         */
    //     }
    // }

    // // Post-processing if needed:
    // // Density and Pressure can be derived here as in MATLAB
}

