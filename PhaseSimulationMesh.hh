#pragma once
#ifndef PHASE_SIMULATION_MESH
#define PHASE_SIMULATION_MESH

#include <iostream>
#include <vector>
#include "PhaseNode.hh"

class PhaseSimulationMesh{
    private:
        std::vector<long int> Dimension{64, 64, 1};
        std::vector<PhaseNode> SimuNode;
    public:
        PhaseSimulationMesh(){};
        PhaseSimulationMesh(PhaseNode Nodes);
        PhaseSimulationMesh(std::vector<long> Dimension, PhaseNode Nodes);
        ~PhaseSimulationMesh(){};

        long getNumNodes(); // Return the Number of Nodes in Mesh
        void fillNodes(PhaseNode Nodes);

        void Ave_Dis_Noise(long double Ave, long double Var);
        PhaseNode& findNode(long X, long Y, long Z);

        void showProperties();
        void showCon();
};

        
PhaseSimulationMesh::PhaseSimulationMesh(PhaseNode Nodes){
    fillNodes(Nodes);
}

PhaseSimulationMesh::PhaseSimulationMesh(std::vector<long> SizeInfo, PhaseNode Nodes){
    Dimension = SizeInfo;
    fillNodes(Nodes);
}

long PhaseSimulationMesh::getNumNodes(){
    return SimuNode.size();
}

void PhaseSimulationMesh::fillNodes(PhaseNode Nodes){
    for(long i = 0; i < Dimension.at(0)*Dimension.at(1)*Dimension.at(2); i++){
        SimuNode.push_back(Nodes);
    }
}

void PhaseSimulationMesh::Ave_Dis_Noise(long double Ave, long double Var){
    for(auto &Nodes : SimuNode){
        Nodes.ConInitial_AveDis(Ave,Var);
    }
}


PhaseNode& PhaseSimulationMesh::findNode(long X, long Y, long Z){
    if(X<Dimension.at(0)&&Y<Dimension.at(1)&&Z<Dimension.at(2) && !(X<0) &&!(Y<0) &&!(Z<0)){
        return SimuNode.at(X+Y*Dimension.at(0)+Z*Dimension.at(0)*Dimension.at(1));
    }
    else throw std::invalid_argument("Index Not in Mesh");
}

void PhaseSimulationMesh::showProperties(){
    std::cout<<"PhaseSimulationMesh Properties:"<<std::endl;
    std::cout<<"Mesh Size:\t\t"<<Dimension.at(0)<<"\u0078"<<Dimension.at(1)<<"\u0078"<<Dimension.at(2)<<std::endl;
    std::cout<<"Number of Nodes:\t"<<getNumNodes()<<std::endl;
}

void PhaseSimulationMesh::showCon(){
    for(long i = 0; i < getNumNodes(); i++){
        SimuNode.at(i).showNode();

    }
}

#endif