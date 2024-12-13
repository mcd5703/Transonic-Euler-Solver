#ifndef HLLCFLUX_HPP
#define HLLCFLUX_HPP

#include <vector>
#include <cmath>

// HLLCFlux function computes the HLLC flux given left and right states (qL, qR),
// the gamma value (gamInf), and the face geometry (xx, xy, vol).
// qL, qR: [rho, rho*u, rho*v, rho*E]
// gamInf: ratio of specific heats (gamma)
// xx, xy: face orientation vectors
// vol: volume scale factor (may be used in normalization)
// F: output flux vector (size 4)
inline void HLLCFlux(const std::vector<double> &qL, const std::vector<double> &qR,
                        double gamInf, double xx, double xy, double vol,
                        std::vector<double> &F)
{
    double xs = std::sqrt(xx*xx + xy*xy);
    double nx = xx/xs;
    double ny = xy/xs;

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
    double u = (std::sqrt(rhoL)*uL + std::sqrt(rhoR)*uR)/(std::sqrt(rhoL) + std::sqrt(rhoR));
    double vv = (std::sqrt(rhoL)*vL + std::sqrt(rhoR)*vR)/(std::sqrt(rhoL) + std::sqrt(rhoR));
    double H = (std::sqrt(rhoL)*HL + std::sqrt(rhoR)*HR)/(std::sqrt(rhoL) + std::sqrt(rhoR));
    double T = (gamInf-1)*(H - 0.5*(u*u + vv*vv))/gamInf;
    double c = std::sqrt(gamInf*T);
    double UU = u*nx + vv*ny;
    double SL = UU - c;
    double SR = UU + c;

    // Compute wave speeds and states
    double R1 = SL - UUL;
    double R2 = SR - UUR;
    double SS = (rhoR*TR - rhoL*TL + rhoL*UUL*(SL - UUL) - rhoR*UUR*(SR - UUR)) /
                (rhoL*R1 - rhoR*R2);

    std::vector<double> DS(4,0.0);
    DS[0]=0; DS[1]=nx; DS[2]=ny; DS[3]=SS;

    if (SL > 0.0) {
        // Supersonic left
        F = FL;
    } else if (SR < 0.0) {
        // Supersonic right
        F = FR;
    } else if (SS > 0.0) {
        // Left intermediate state
        std::vector<double> tmp(4,0.0);
        for (int n=0; n<4; n++){
            double A = (SL*qL[n]-FL[n])*SS;
            double B = SL*(rhoL*TL + rhoL*(SL - UUL)*(SS - UUL))*DS[n];
            tmp[n] = (A + B)/(SL - SS);
        }
        F = tmp;
    } else {
        // Right intermediate state
        std::vector<double> tmp(4,0.0);
        for (int n=0; n<4; n++){
            double A = (SR*qR[n]-FR[n])*SS;
            double B = SR*(rhoR*TR + rhoR*(SR - UUR)*(SS - UUR))*DS[n];
            tmp[n] = (A + B)/(SR - SS);
        }
        F = tmp;
    }

    // Multiply by face area
    for (int n=0; n<4; n++)
        F[n] *= xs;
}

#endif // HLLCFLUX_HPP
