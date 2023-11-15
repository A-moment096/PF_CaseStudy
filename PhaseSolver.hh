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
    std::vector<double> Laplacian (int whichSTNCL, PhaseSimulationMesh &mesh, int whichPara){
        std::vector<double> result(mesh.getNumNodes(),0.0);
        double dx = mesh.getStepLength(WHICHDIM::X);
        double dy = mesh.getStepLength(WHICHDIM::Y);
        double dz = mesh.getStepLength(WHICHDIM::Z);
        if(whichSTNCL = STENCILE::FIVEPOINT){
#pragma omp parallel for
            for(int i = 0; i < mesh.getNumNodes(); i++){
                
                double c = mesh.findNode(i).getProperties(whichPara).at(0);
                double f = mesh.findNode(i).Forward->getProperties(whichPara).at(0);
                double b = mesh.findNode(i).Backward->getProperties(whichPara).at(0);
                double l = mesh.findNode(i).Left->getProperties(whichPara).at(0);
                double r = mesh.findNode(i).Right->getProperties(whichPara).at(0);

                result.at(i) = ((f+b+l+r-4*c)/(dx*dy*dz));
            }
        }
        return result;
    }

};


#endif