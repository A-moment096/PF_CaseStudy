#pragma once
#ifndef SPINODAL_SOLVER
#define SPINODAL_SOLVER

#include <iostream>
#include <vector>
#include "PhaseDomain.hh"

class SpinodalSolver: public PhaseDomain
{
private:
    /* data */
public:
    SpinodalSolver(/* args */);
    ~SpinodalSolver();
};

SpinodalSolver::SpinodalSolver(/* args */)
{
}

SpinodalSolver::~SpinodalSolver()
{
}


#endif