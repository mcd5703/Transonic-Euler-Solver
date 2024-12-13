#include "solver.hpp"
#include "grid.hpp"
#include "boundaryconditions.hpp"
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

}