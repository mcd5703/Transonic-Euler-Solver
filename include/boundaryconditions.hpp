#ifndef BOUNDARYCONDITIONS_HPP
#define BOUNDARYCONDITIONS_HPP

#include <vector>

// bcwall returns a new qB array after applying the wall boundary condition.
std::vector<std::vector<std::vector<double>>> bcwall(
    const std::vector<std::vector<std::vector<double>>> &q,
    const std::vector<double> &qInfPrim,
    const std::vector<std::vector<double>> &xx,
    const std::vector<std::vector<double>> &xy,
    const std::vector<std::vector<double>> &vol,
    int NJ, int NK);

// bcfree returns a new qB array after applying the free-stream boundary condition.
std::vector<std::vector<std::vector<double>>> bcfree(
    const std::vector<std::vector<std::vector<double>>> &q,
    const std::vector<double> &qInfPrim,
    const std::vector<std::vector<double>> &xx,
    const std::vector<std::vector<double>> &xy,
    const std::vector<std::vector<double>> &vol,
    int NJ, int NK);

#endif // BOUNDARYCONDITIONS_HPP
