#pragma once
#ifndef PHASE_DOMAIN
#define PHASE_DOMAIN

#include <iostream>
#include <vector>
#include "PhaseNode.hh"
#include "PhaseSimulationMesh.hh"

class PhaseDomain{
    private:
        std::vector<PhaseNode> GrainsNodes{DefaultNode}; // nodes for different grains
        PhaseSimulationMesh BaseMesh = DefaultMesh; // whole base mesh
    public:
        PhaseDomain(){};
        PhaseDomain(std::vector<PhaseNode> SimuNodes, PhaseSimulationMesh SimuGrains); // node List and base mesh
        ~PhaseDomain(){};

        // PhaseNode& findNode(long X, long Y, long Z);

        void showGrainProp();
        void showBaseMeshProp();
} DefaultDomain ;

PhaseDomain::PhaseDomain(std::vector<PhaseNode> SimuNodes, PhaseSimulationMesh SimuGrains){
   GrainsNodes = SimuNodes;
   BaseMesh = SimuGrains;
}

void PhaseDomain::showGrainProp(){
    std::cout<<"Grain Number:\t"<<GrainsNodes.size()<<std::endl;
    std::cout<<"************************************************"<<std::endl;
    int order = 0;
    for(auto grain : GrainsNodes){
        std::cout<<"Grain No. "<<order<<std::endl;
        grain.showNode();
        std::cout<<std::endl;
        order++;
    }
}

void PhaseDomain::showBaseMeshProp(){
    BaseMesh.showMeshProp();
}

#endif