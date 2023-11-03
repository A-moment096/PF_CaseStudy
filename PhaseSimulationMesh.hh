#pragma once
#ifndef PHASE_SIMULATION_MESH
#define PHASE_SIMULATION_MESH

#include <iostream>
#include <vector>
#include "PhaseNode.hh"

class PhaseSimulationMesh{
    private:
        std::vector<long int> Dimension{64, 64, 1};
        long int& BoxX = Dimension.at(0);
        long int& BoxY = Dimension.at(1);
        long int& BoxZ = Dimension.at(2);

        std::vector<PhaseNode> SimuNodes;

    public:
        PhaseSimulationMesh(){};
        PhaseSimulationMesh(PhaseNode Nodes);
        PhaseSimulationMesh(std::vector<long> SizeInfo, PhaseNode Nodes);
        ~PhaseSimulationMesh(){};

        long getNumNodes(); // Return the Number of Nodes in Mesh
        void fillNodes(PhaseNode Nodes);

        void Ave_Dis_Noise(long double Ave, long double Var);
        PhaseNode& findNode(long X, long Y, long Z);

        void showMeshProp();
        void showNodeProp(unsigned index);
};

        
PhaseSimulationMesh::PhaseSimulationMesh(PhaseNode Nodes){
    fillNodes(Nodes);
}

PhaseSimulationMesh::PhaseSimulationMesh(std::vector<long> SizeInfo, PhaseNode Nodes){
    Dimension = SizeInfo;
    fillNodes(Nodes);
}

long PhaseSimulationMesh::getNumNodes(){
    return SimuNodes.size();
}

void PhaseSimulationMesh::fillNodes(PhaseNode Nodes){
    for(long i = 0; i < BoxX*BoxY*BoxZ; i++){
        SimuNodes.push_back(Nodes);
    }
}

void PhaseSimulationMesh::Ave_Dis_Noise(long double Ave, long double Var){
    for(auto &Nodes : SimuNodes){
        Nodes.ConInitial_AveDis(Ave,Var);
    }
}


PhaseNode& PhaseSimulationMesh::findNode(long X, long Y, long Z){
    if(X<BoxX&&Y<BoxY&&Z<BoxY && !(X<0) &&!(Y<0) &&!(Z<0)){
        return SimuNodes.at(X+Y*BoxX+Z*BoxX*BoxY);
    }
    else throw std::invalid_argument("Index Not in Mesh");
}

void PhaseSimulationMesh::showMeshProp(){
    std::cout<<"PhaseSimulationMesh Properties:"<<std::endl;
    std::cout<<"Mesh Size:\t\t"<<BoxX<<"\u0078"<<BoxY<<"\u0078"<<BoxY<<std::endl;
    std::cout<<"Number of Nodes:\t"<<getNumNodes()<<std::endl;
}

void PhaseSimulationMesh::showNodeProp(unsigned index){
    for(long i = 0; i < getNumNodes(); i++){
        std::cout<<SimuNodes.at(i).getProperties(index);
        if(i%BoxX==BoxX-1)std::cout<<std::endl;
        if(i%BoxX*BoxY == BoxX*BoxY-1)std::cout<<"\n"<<std::endl;
    }
}

#endif