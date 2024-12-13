#include "HLLCFlux.hpp"
#include <cmath>
#include <vector>

void HLLCFlux(
    const std::vector<double> &qL,
    const std::vector<double> &qR,
    double gamInf,
    double xx,
    double xy,
    double vol,
    std::vector<double> &F
)
{
    double xs = std::sqrt(xx*xx + xy*xy); // Face area
    double nx = xx/xs;
    double ny = xy/xs;

    // Left state
    double rhoL = qL[0];
    double uL = qL[1]/rhoL;
    double vL = qL[2]/rhoL;
    double TL = (gamInf-1)*(qL[3]/rhoL - 0.5*(uL*uL + vL*vL));
    double HL = gamInf*TL/(gamInf-1) + 0.5*(uL*uL + vL*vL);
    double UUL = uL*nx + vL*ny;

    std::vector<double> FL(4,0.0);
    FL[0] = rhoL*UUL;
    FL[1] = rhoL*UUL*uL + rhoL*TL*nx;
    FL[2] = rhoL*UUL*vL + rhoL*TL*ny;
    FL[3] = rhoL*UUL*HL;

    // Right state
    double rhoR = qR[0];
    double uR = qR[1]/rhoR;
    double vR = qR[2]/rhoR;
    double TR = (gamInf-1)*(qR[3]/rhoR - 0.5*(uR*uR + vR*vR));
    double HR = gamInf*TR/(gamInf-1) + 0.5*(uR*uR + vR*vR);
    double UUR = uR*nx + vR*ny;

    std::vector<double> FR(4,0.0);
    FR[0] = rhoR*UUR;
    FR[1] = rhoR*UUR*uR + rhoR*TR*nx;
    FR[2] = rhoR*UUR*vR + rhoR*TR*ny;
    FR[3] = rhoR*UUR*HR;

    double rho = std::sqrt(rhoL*rhoR);
    double u = (std::sqrt(rhoL)*uL + std::sqrt(rhoR)*uR)/(std::sqrt(rhoL)+std::sqrt(rhoR));
    double vv = (std::sqrt(rhoL)*vL + std::sqrt(rhoR)*vR)/(std::sqrt(rhoL)+std::sqrt(rhoR));
    double H = (std::sqrt(rhoL)*HL + std::sqrt(rhoR)*HR)/(std::sqrt(rhoL)+std::sqrt(rhoR));
    // T used to define speed of sound c
    // T = (gamInf-1)*(H - 0.5*(u^2+v^2))/gamInf; but we only need c:
    double c = std::sqrt((gamInf-1)*(H - 0.5*(u*u + vv*vv)));

    double UU = u*nx + vv*ny;
    double SL = UU - c;
    double SR = UU + c;

    // Compute SS (contact wave speed)
    double numerator = (rhoR*TR - rhoL*TL + rhoL*UUL*(SL - UUL) - rhoR*UUR*(SR - UUR));
    double denominator = (rhoL*(SL - UUL) - rhoR*(SR - UUR));
    double SS = numerator/denominator;

    double PLR = 0.5*(rhoL*TL + rhoR*TR + rhoL*(SL - UUL)*(SS - UUL) + rhoR*(SR - UUR)*(SS - UUR));
    std::vector<double> DS(4,0.0);
    DS[0] = 0; DS[1] = nx; DS[2] = ny; DS[3] = SS;

    std::vector<double> Ftmp(4,0.0);

    if (SL > 0) {
        // Supersonic flow from left to right
        Ftmp = FL;
    } else if (SR < 0) {
        // Supersonic flow from right to left
        Ftmp = FR;
    } else if (SS > 0) {
        // Left intermediate state
        for (int n=0; n<4; n++) {
            Ftmp[n] = (SS*(SL*qL[n] - FL[n]) + SL*(rhoL*TL + rhoL*(SL - UUL)*(SS - UUL))*DS[n])/(SL - SS);
        }
    } else {
        // Right intermediate state
        for (int n=0; n<4; n++) {
            Ftmp[n] = (SS*(SR*qR[n]-FR[n]) + SR*(rhoR*TR + rhoR*(SR - UUR)*(SS - UUR))*DS[n])/(SR - SS);
        }
    }

    // Multiply by face area
    for (int n=0; n<4; n++){
        Ftmp[n] = xs * Ftmp[n];
    }

    F = Ftmp;
}
