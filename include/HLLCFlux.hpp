#ifndef HLLCFLUX_HPP
#define HLLCFLUX_HPP

#include <vector>

// HLLCFlux computes the HLLC flux given left and right states qL and qR,
// gamma (gamInf), and geometric parameters xx, xy, and vol.
// qL and qR are size 4: [rho, rho*u, rho*v, rho*E].
// Output F is also size 4.
void HLLCFlux(
    const std::vector<double> &qL,
    const std::vector<double> &qR,
    double gamInf,
    double xx,
    double xy,
    double vol,
    std::vector<double> &F
);

#endif // HLLCFLUX_HPP
