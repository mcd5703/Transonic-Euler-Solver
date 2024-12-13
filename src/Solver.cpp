#include "solver.hpp"
#include "grid.hpp"
#include "boundaryconditions.hpp"
#include "simpleFunctions.hpp"
#include "resid.hpp"
#include "lusgs.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include <matplot/matplot.h>

void visualizeGrid(const std::vector<std::vector<double>> &x,
                    const std::vector<std::vector<double>> &y,
                    int NJ, int NK)
{
    using namespace matplot;

    // Visualization
    auto f = figure(true);
    f->name("Grid Plot");
    f->number_title(false);
    f->color("black");
    f->position({0, 0, 600, 600});
    f->size(1000, 1000);
    f->font("Arial");
    f->font_size(30);
    f->title("Grid Plot");

    // Get current axes
    auto ax = f->current_axes();
    ax->hold(true);

    // Plot lines in J direction (constant j)
    for (int j=0;j<NJ;j++){
        std::vector<double> xv, yv;
        xv.reserve(NK);
        yv.reserve(NK);
        for (int k=0;k<NK;k++){
            xv.push_back(x[j][k]);
            yv.push_back(y[j][k]);
        }
        ax->plot(xv,yv,"b");
    }

    // Plot lines in K direction (constant k)
    for (int k=0;k<NK;k++){
        std::vector<double> xv, yv;
        xv.reserve(NJ);
        yv.reserve(NJ);
        for (int j=0;j<NJ;j++){
            xv.push_back(x[j][k]);
            yv.push_back(y[j][k]);
        }
        ax->plot(xv,yv,"r");
    }

    // Set axis limits and labels
    ax->xlim({-1,2});
    ax->ylim({-1,1});

    show();
}

void Solver(int NACA, double Minf, double alpha, double gamInf,
            int itmax, double CFL, int nDisp, int NJ, int NK, double R2)
{

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

    // Generate the grid
    std::vector<std::vector<double>> x(NJ, std::vector<double>(NK,0.0));
    std::vector<std::vector<double>> y(NJ, std::vector<double>(NK,0.0));
    makegrid(NACA,R2,NJ,NK,x,y);

    // Visualize the grid before proceeding
    visualizeGrid(x,y,NJ,NK);

    // Compute metrics
    std::vector<std::vector<double>> xx(NJ, std::vector<double>(NK,0.0));
    std::vector<std::vector<double>> xy(NJ, std::vector<double>(NK,0.0));
    std::vector<std::vector<double>> yx(NJ, std::vector<double>(NK,0.0));
    std::vector<std::vector<double>> yy(NJ, std::vector<double>(NK,0.0));
    std::vector<std::vector<double>> vol(NJ, std::vector<double>(NK,0.0));
    metrics(x,y,NJ,NK,xx,xy,yx,yy,vol);

    // Initialize solution to free stream
    double rhoInf = 1.0; // Density Normalized to 1
    double uInf = Minf * std::cos(alpha * M_PI/180.0); 
    double vInf = Minf * std::sin(alpha * M_PI/180.0);
    double TInf = 1.0/gamInf; // Temperature normalized by gamma
    std::vector<double> qInfPrim = {rhoInf, uInf, vInf, TInf, gamInf};

    std::vector<std::vector<std::vector<double>>> q(NJ,
        std::vector<std::vector<double>>(NK,std::vector<double>(4,0.0)));

    for (int j = 0; j < NJ; j++) {
        for (int k = 0; k < NK; k++) {
            q[j][k][0] = rhoInf; // Density
            q[j][k][1] = rhoInf*uInf; // XMomentum
            q[j][k][2] = rhoInf*vInf; // YMomentum
            q[j][k][3] = rhoInf*(TInf/(gamInf-1.0) + 0.5*(uInf*uInf + vInf*vInf));
        }
    }

    std::vector<std::vector<std::vector<double>>> qn = q;
    std::vector<std::vector<std::vector<double>>> R(NJ,std::vector<std::vector<double>>(NK,std::vector<double>(4,0.0)));
    std::vector<std::vector<double>> dt(NJ,std::vector<double>(NK,0.0));

    // We'll need v1, v2 arrays same dimension as q:
    std::vector<std::vector<std::vector<double>>> v1(NJ,std::vector<std::vector<double>>(NK,std::vector<double>(4,0.0)));
    std::vector<std::vector<std::vector<double>>> v2(NJ,std::vector<std::vector<double>>(NK,std::vector<double>(4,0.0)));
    std::vector<std::vector<std::vector<double>>> RN(NJ,std::vector<std::vector<double>>(NK,std::vector<double>(4,0.0)));

    // Arrays for Mplot, P, Cp
    std::vector<std::vector<double>> Mplot(NJ,std::vector<double>(NK,0.0));
    std::vector<std::vector<double>> Pmat(NJ,std::vector<double>(NK,0.0));
    std::vector<std::vector<double>> Cp(NJ,std::vector<double>(NK,0.0));

    using namespace matplot;
    auto f = figure(true);
    f->name("Flow Field");
    f->number_title(false);
    f->color("black");
    f->position({0, 0, 600, 600});
    f->size(1000, 1000);
    f->font("Arial");
    f->font_size(30);
    f->title("Mach Contours Evolution");

    auto ax = f->current_axes();
    ax->hold(true);

    for (int it=1; it<=itmax; it++) {
        // q(1,:,:) = 1.5*q(2,:,:) - 0.5*q(3,:,:)
        for (int k=0; k<NK; k++){
            for (int n=0; n<4; n++){
                q[0][k][n] = 1.5*q[1][k][n] - 0.5*q[2][k][n];
            }
        }

        // bcwall on top boundary (j=0)
        {
            std::vector<std::vector<std::vector<double>>> qSlice(1,std::vector<std::vector<double>>(NK,std::vector<double>(4,0.0)));
            for (int k=0; k<NK; k++){
                for (int n=0; n<4; n++){
                    qSlice[0][k][n] = q[0][k][n];
                }
            }

            // Make single-line xx,xy,vol arrays
            std::vector<std::vector<double>> xxLine(1,std::vector<double>(NK));
            std::vector<std::vector<double>> xyLine(1,std::vector<double>(NK));
            std::vector<std::vector<double>> volLine(1,std::vector<double>(NK));
            for (int kk=0; kk<NK; kk++){
                xxLine[0][kk]=xx[0][kk];
                xyLine[0][kk]=xy[0][kk];
                volLine[0][kk]=vol[0][kk];
            }

            auto qB = bcwall(qSlice, qInfPrim, xxLine, xyLine, volLine, 1, NK);
            for (int k=0; k<NK; k++){
                for (int n=0; n<4; n++){
                    q[0][k][n] = qB[0][k][n];
                }
            }
        }

        // q(NJ,:,:) = q(NJ-1,:,:)
        for (int k=0; k<NK; k++){
            for (int n=0; n<4; n++){
                q[NJ-1][k][n] = q[NJ-2][k][n];
            }
        }

        // bcfree on bottom boundary (j=NJ-1)
        {
            std::vector<std::vector<std::vector<double>>> qSliceBottom(1,std::vector<std::vector<double>> (NK,std::vector<double>(4,0.0)));
            for (int k=0;k<NK;k++){
                for (int n=0;n<4;n++){
                    qSliceBottom[0][k][n]=q[NJ-1][k][n];
                }
            }

            std::vector<std::vector<double>> xxBottomLine(1,std::vector<double>(NK));
            std::vector<std::vector<double>> xyBottomLine(1,std::vector<double>(NK));
            std::vector<std::vector<double>> volBottomLine(1,std::vector<double>(NK));

            for (int kk = 0; kk < NK; kk++){
                xxBottomLine[0][kk] = xx[NJ-1][kk];
                xyBottomLine[0][kk] = xy[NJ-1][kk];
                volBottomLine[0][kk] = vol[NJ-1][kk];
            }

            auto qBbottom = bcfree(qSliceBottom, qInfPrim, xxBottomLine, xyBottomLine, volBottomLine, 1, NK);
            for (int k=0;k<NK;k++){
                for (int n=0;n<4;n++){
                    q[NJ-1][k][n]=qBbottom[0][k][n];
                }
            }
        }

        qn = q;

        resid(q,xx,xy,yx,yy,vol,qInfPrim,NJ,NK,CFL,R,dt);
        lusgs(R,q,xx,xy,yx,yy,vol,dt,gamInf,NJ,NK,v1);
        RN = R;

        double h11 = qdotq(v1,v1,NJ,NK);
        h11 = std::sqrt(h11);
        for (int j=0;j<NJ;j++){
            for (int k=0;k<NK;k++){
                for (int n=0;n<4;n++){
                    v1[j][k][n] /= h11;
                }
            }
        }

        double eps = 2e-8;
        for (int j=0;j<NJ;j++){
            for (int k=0;k<NK;k++){
                for (int n=0;n<4;n++){
                    q[j][k][n] = qn[j][k][n] + eps*v1[j][k][n];
                }
            }
        }

        // q(1,:,:) = 1.5*q(2,:,:) - 0.5*q(3,:,:)
        for (int k=0;k<NK;k++){
            for (int n=0;n<4;n++){
                q[0][k][n] = 1.5*q[1][k][n] - 0.5*q[2][k][n];
            }
        }

        // bcwall top boundary again
        {
            std::vector<std::vector<std::vector<double>>> qSlice(1,std::vector<std::vector<double>> (NK,std::vector<double>(4,0.0)));
            for (int kk=0; kk<NK; kk++){
                for (int n=0;n<4;n++){
                    qSlice[0][kk][n]=q[0][kk][n];
                }
            }

            std::vector<std::vector<double>> xxLine(1,std::vector<double>(NK));
            std::vector<std::vector<double>> xyLine(1,std::vector<double>(NK));
            std::vector<std::vector<double>> volLine(1,std::vector<double>(NK));
            for (int kk=0; kk<NK; kk++){
                xxLine[0][kk]=xx[0][kk];
                xyLine[0][kk]=xy[0][kk];
                volLine[0][kk]=vol[0][kk];
            }

            auto qB=bcwall(qSlice,qInfPrim,xxLine,xyLine,volLine,1,NK);
            for (int kk=0;kk<NK;kk++){
                for (int n=0;n<4;n++){
                    q[0][kk][n]=qB[0][kk][n];
                }
            }
        }

        // q(NJ,:,:) = q(NJ-1,:,:)
        for (int k=0;k<NK;k++){
            for (int n=0;n<4;n++){
                q[NJ-1][k][n]=q[NJ-2][k][n];
            }
        }

        // bcfree bottom boundary again
        {
            std::vector<std::vector<std::vector<double>>> qSliceBottom(1,std::vector<std::vector<double>> (NK,std::vector<double>(4,0.0)));
            for (int k=0;k<NK;k++){
                for (int n=0;n<4;n++){
                    qSliceBottom[0][k][n]=q[NJ-1][k][n];
                }
            }

            std::vector<std::vector<double>> xxBottomLine(1,std::vector<double>(NK));
            std::vector<std::vector<double>> xyBottomLine(1,std::vector<double>(NK));
            std::vector<std::vector<double>> volBottomLine(1,std::vector<double>(NK));
            for (int kk=0; kk<NK; kk++){
                xxBottomLine[0][kk]=xx[NJ-1][kk];
                xyBottomLine[0][kk]=xy[NJ-1][kk];
                volBottomLine[0][kk]=vol[NJ-1][kk];
            }

            auto qB=bcfree(qSliceBottom,qInfPrim,xxBottomLine,xyBottomLine,volBottomLine,1,NK);
            for (int k=0;k<NK;k++){
                for (int n=0;n<4;n++){
                    q[NJ-1][k][n]=qB[0][k][n];
                }
            }
        }

        resid(q,xx,xy,yx,yy,vol,qInfPrim,NJ,NK,CFL,R,dt);
        for (int j=0;j<NJ;j++){
            for (int k=0;k<NK;k++){
                for (int n=0;n<4;n++){
                    v2[j][k][n] = v1[j][k][n] + (R[j][k][n] - RN[j][k][n])/eps;
                }
            }
        }

        double h22 = qdotq(v2,v2,NJ,NK);
        double h12 = qdotq(R,v2,NJ,NK);
        double a1 = -h12/h22;
        a1 = std::max(a1,0.5*h11);

        for (int j=0;j<NJ;j++){
            for (int k=0;k<NK;k++){
                for (int n=0;n<4;n++){
                    q[j][k][n] = qn[j][k][n] + a1*v1[j][k][n];
                }
            }
        }

        if ((it % nDisp) == 0) {
            double sumsq=0.0;
            int count=0;
            for (int jj=1;jj<NJ-1;jj++){
                for (int kk=0;kk<NK;kk++){
                    sumsq += R[jj][kk][0]*R[jj][kk][0];
                    count++;
                }
            }
            double rhoresid = std::sqrt(sumsq)/CFL;

            for (int j=0;j<NJ;j++){
                for (int k=0;k<NK;k++){
                    double rho=q[j][k][0];
                    double rhou=q[j][k][1];
                    double rhov=q[j][k][2];
                    double rhoE=q[j][k][3];
                    double u=rhou/rho;
                    double v=rhov/rho;
                    double p=(gamInf-1.0)*(rhoE -0.5*(rhou*rhou+rhov*rhov)/rho);
                    double a2 = gamInf*(gamInf-1)*(rhoE/rho -0.5*(u*u+v*v));
                    double M = std::sqrt((u*u+v*v)/a2);
                    Mplot[j][k]=M;
                    Pmat[j][k]=p;
                    Cp[j][k]=(p -1.0/gamInf)/(0.5*Minf*Minf);
                }
            }

            double cx=0.0, cy=0.0;
            int JJ=0;
            for (int K=0; K<NK-1; K++){
                cx -= 0.5*(Cp[JJ][K]*xx[JJ][K] + Cp[JJ][K+1]*xx[JJ][K+1]);
                cy -= 0.5*(Cp[JJ][K]*xy[JJ][K] + Cp[JJ][K+1]*xy[JJ][K+1]);
            }

            double cl = cy*std::cos(alpha*M_PI/180.0) - cx*std::sin(alpha*M_PI/180.0);
            double cd = cx*std::cos(alpha*M_PI/180.0) + cy*std::sin(alpha*M_PI/180.0);

            std::cout << "Iter " << it << ": resid = " << rhoresid << ", cl = " << cl << ", cd = " << cd << std::endl;

            ax->clear();
            ax->contourf(x, y, Mplot, 21);
            ax->xlim({-1,2});
            ax->ylim({-1,1});
            f->draw();
        }
    }

    // Density = q[:,:,1]*0.97 (C++: just scale q[j][k][0]*=0.97 if desired)
    // Pressure = P*70989.45*1.4; (You can compute final Pressure similarly if needed)

}
