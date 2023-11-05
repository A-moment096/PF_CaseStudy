#pragma once
#ifndef PHASE_GRAIN
#define PHASE_GRAIN

#include <iostream>
#include <vector>
#include "PhaseNode.hh"
#include "PhaseSimulationMesh.hh"

class PhaseGrain{
    private:
        unsigned NumGrains = 1;
        std::vector<PhaseSimulationMesh> Grains;
        //long int& BoxX = Dimension.at(0);
        //long int& BoxY = Dimension.at(1);
        //long int& BoxZ = Dimension.at(2);

    public:
        PhaseGrain(){};
        PhaseGrain(PhaseSimulationMesh mesh);
        PhaseGrain(unsigned num, PhaseSimulationMesh mesh);
        ~PhaseGrain(){};

        
};

PhaseGrain::PhaseGrain(PhaseSimulationMesh mesh){
    Grains.push_back(mesh);
}    

PhaseGrain::PhaseGrain(unsigned num, PhaseSimulationMesh mesh){
    for(unsigned i = 0; i < num; i++){
        Grains.push_back(mesh);
    }
}

#endif