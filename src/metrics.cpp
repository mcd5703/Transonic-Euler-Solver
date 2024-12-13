#include "grid.hpp"
#include <vector>

void metrics(const std::vector<std::vector<double>> &x,
                const std::vector<std::vector<double>> &y,
                int NJ, int NK,
                std::vector<std::vector<double>> &xx,
                std::vector<std::vector<double>> &xy,
                std::vector<std::vector<double>> &yx,
                std::vector<std::vector<double>> &yy,
                std::vector<std::vector<double>> &vol)
{
    std::vector<std::vector<double>> xxsi(NJ, std::vector<double>(NK,0.0));
    std::vector<std::vector<double>> yxsi(NJ, std::vector<double>(NK,0.0));
    std::vector<std::vector<double>> xeta(NJ, std::vector<double>(NK,0.0));
    std::vector<std::vector<double>> yeta(NJ, std::vector<double>(NK,0.0));

    // Compute xxsi, yxsi
    for (int k=0; k<NK; k++) {
        xxsi[0][k] = x[1][k] - x[0][k];
        yxsi[0][k] = y[1][k] - y[0][k];
        for (int j=1; j<NJ-1; j++) {
            xxsi[j][k] = 0.5*(x[j+1][k]-x[j-1][k]);
            yxsi[j][k] = 0.5*(y[j+1][k]-y[j-1][k]);
        }
        xxsi[NJ-1][k] = x[NJ-1][k] - x[NJ-2][k];
        yxsi[NJ-1][k] = y[NJ-1][k] - y[NJ-2][k];
    }

    // Compute xeta, yeta
    for (int j=0; j<NJ; j++) {
        xeta[j][0] = 0.5*(x[j][1] - x[j][NK-2]);
        yeta[j][0] = 0.5*(y[j][1] - y[j][NK-2]);
        for (int k=1; k<NK-1; k++) {
            xeta[j][k] = 0.5*(x[j][k+1]-x[j][k-1]);
            yeta[j][k] = 0.5*(y[j][k+1]-y[j][k-1]);
        }
        xeta[j][NK-1] = xeta[j][0];
        yeta[j][NK-1] = yeta[j][0];
    }

    // Compute vol, xx, xy, yx, yy
    for (int j=0; j<NJ; j++) {
        for (int k=0; k<NK; k++) {
            vol[j][k] = xxsi[j][k]*yeta[j][k] - xeta[j][k]*yxsi[j][k];
            xx[j][k] = yeta[j][k];
            xy[j][k] = -xeta[j][k];
            yx[j][k] = -yxsi[j][k];
            yy[j][k] = xxsi[j][k];
        }
    }
}
