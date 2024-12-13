#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <algorithm>

// Use Eigen for vectors and matrices
using namespace Eigen;
typedef Vector4d Vector4;

// Function: HLLCFlux
Vector4 HLLCFlux(const Vector4 &qL, const Vector4 &qR, double gamInf, double xx, double xy, double vol) {
    // Compute face normal and area
    double xs = std::sqrt(xx * xx + xy * xy); // Face area
    double nx = xx / xs;                     // Normalized x-component
    double ny = xy / xs;                     // Normalized y-component

    // Left state variables
    double rhoL = qL(0);
    double uL = qL(1) / rhoL;
    double vL = qL(2) / rhoL;
    double TL = (gamInf - 1.0) * (qL(3) / rhoL - 0.5 * (uL * uL + vL * vL));
    double HL = gamInf * TL / (gamInf - 1.0) + 0.5 * (uL * uL + vL * vL);
    double UUL = uL * nx + vL * ny; // Normal velocity
    Vector4 FL;
    FL << rhoL * UUL,
          rhoL * UUL * uL + rhoL * TL * nx,
          rhoL * UUL * vL + rhoL * TL * ny,
          rhoL * UUL * HL;

    // Right state variables
    double rhoR = qR(0);
    double uR = qR(1) / rhoR;
    double vR = qR(2) / rhoR;
    double TR = (gamInf - 1.0) * (qR(3) / rhoR - 0.5 * (uR * uR + vR * vR));
    double HR = gamInf * TR / (gamInf - 1.0) + 0.5 * (uR * uR + vR * vR);
    double UUR = uR * nx + vR * ny; // Normal velocity
    Vector4 FR;
    FR << rhoR * UUR,
          rhoR * UUR * uR + rhoR * TR * nx,
          rhoR * UUR * vR + rhoR * TR * ny,
          rhoR * UUR * HR;

    // Roe-averaged quantities
    double sqrtRhoL = std::sqrt(rhoL);
    double sqrtRhoR = std::sqrt(rhoR);
    double rho = sqrtRhoL * sqrtRhoR;
    double u = (sqrtRhoL * uL + sqrtRhoR * uR) / (sqrtRhoL + sqrtRhoR);
    double v = (sqrtRhoL * vL + sqrtRhoR * vR) / (sqrtRhoL + sqrtRhoR);
    double H = (sqrtRhoL * HL + sqrtRhoR * HR) / (sqrtRhoL + sqrtRhoR);
    double c = std::sqrt((gamInf - 1.0) * (H - 0.5 * (u * u + v * v)));

    // Wave speeds
    double UU = u * nx + v * ny; // Roe-averaged normal velocity
    double SL = UU - c;          // Left wave speed
    double SR = UU + c;          // Right wave speed

    // Intermediate wave speed
    double SS = (rhoR * TR - rhoL * TL + rhoL * UUL * (SL - UUL) - rhoR * UUR * (SR - UUR)) /
                (rhoL * (SL - UUL) - rhoR * (SR - UUR));

    // Intermediate pressure
    double PLR = 0.5 * (rhoL * TL + rhoR * TR + rhoL * (SL - UUL) * (SS - UUL) + rhoR * (SR - UUR) * (SS - UUR));

    // Resultant flux vector
    Vector4 F;

    // HLLC Flux Cases
    if (SL > 0) { // Supersonic flow from left to right
        F = FL;
    } else if (SR < 0) { // Supersonic flow from right to left
        F = FR;
    } else if (SS > 0) { // Left intermediate state
        Vector4 DS;
        DS << 0.0, nx, ny, SS;
        for (int n = 0; n < 4; ++n) {
            F(n) = (SS * (SL * qL(n) - FL(n)) + SL * (rhoL * TL + rhoL * (SL - UUL) * (SS - UUL)) * DS(n)) / (SL - SS);
        }
    } else { // Right intermediate state
        Vector4 DS;
        DS << 0.0, nx, ny, SS;
        for (int n = 0; n < 4; ++n) {
            F(n) = (SS * (SR * qR(n) - FR(n)) + SR * (rhoR * TR + rhoR * (SR - UUR) * (SS - UUR)) * DS(n)) / (SR - SS);
        }
    }

    // Multiply flux by face area
    return xs * F;
}
