#pragma once
#ifndef PHASE_SOLVER
#define PHASE_SOLVER

#include <iostream>
#include <vector>
#include "PhaseNode.hh"
#include "PhaseSimulationMesh.hh"
#include "PhaseGrain.hh"

class PhaseSolver: public PhaseGrain
{
private:
    /* data */
public:
    PhaseSolver();
    ~PhaseSolver();
};

PhaseSolver::PhaseSolver()
{
}

PhaseSolver::~PhaseSolver()
{
}





#endif