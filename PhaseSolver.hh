#pragma once
#ifndef PHASE_SOLVER
#define PHASE_SOLVER

#include <iostream>
#include <vector>
#include "PhaseSimulationMesh.hh"

enum STENCILE{FIVEPOINT=5,NINEPOINT=9};
class PhaseSolver : private PhaseSimulationMesh
{
private:
    /* data */
public:
    PhaseSolver(){}
    ~PhaseSolver(){}
    std::vector<double> Laplacian (int whichSTNCL, PhaseSimulationMesh mesh){
        std::vector<double> result;
        if(whichSTNCL = STENCILE::FIVEPOINT){
            for(int i = 0; i < mesh.getNumNodes(); i++){
                mesh.findNode(i)
            }
        }
    }

};


#endif