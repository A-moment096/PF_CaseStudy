#pragma once
#ifndef PHASE_SIMULATION_MESH
#define PHASE_SIMULATION_MESH

#include <iostream>
#include <vector>
#include "PhaseNode.hh"

class PhaseSimulationMesh{
    private:
        std::vector<int> Dimension{64, 64, 1};
        int& BoxX = Dimension.at(0);
        int& BoxY = Dimension.at(1);
        int& BoxZ = Dimension.at(2);

        std::vector<PhaseNode> SimuNodes;
    public:
        PhaseSimulationMesh(){}; // initial with default size, without nodes
        PhaseSimulationMesh(PhaseNode Nodes); // initial with nodes and default size
        PhaseSimulationMesh(std::vector<int> SizeInfo, PhaseNode Nodes); // initial with size and nodes
        ~PhaseSimulationMesh(){};

        int getNumNodes(); // return the number of nodes in mesh
        void fillNodes(PhaseNode Nodes); // fill mesh with nodes 

/**/    void Ave_Dis_Noise(double Ave, double Var); // add noise to the nodes (concentration) with average distribution
        PhaseNode& findNode(int X, int Y, int Z); // find node in the mesh according to the coordinates

        void showMeshProp(); // show the basic information of the mesh
        void showNodesProp(unsigned index); // show one of the properties of the nodes inside the mesh
};

        
PhaseSimulationMesh::PhaseSimulationMesh(PhaseNode Nodes){
    fillNodes(Nodes);
}

PhaseSimulationMesh::PhaseSimulationMesh(std::vector<int> SizeInfo, PhaseNode Nodes){
    Dimension = SizeInfo;
    fillNodes(Nodes);
}

int PhaseSimulationMesh::getNumNodes(){
    return SimuNodes.size();
}

void PhaseSimulationMesh::fillNodes(PhaseNode Nodes){
    for(int i = 0; i < Dimension.at(0)*Dimension.at(1)*Dimension.at(2); i++){
        SimuNodes.push_back(Nodes);
    }
}

// void PhaseSimulationMesh::Ave_Dis_Noise(double Ave, double Var){
//     for(auto &Nodes : SimuNodes){
//         Nodes.ConInitial_AveDis(Ave,Var);
//     }
// }


PhaseNode& PhaseSimulationMesh::findNode(int X, int Y, int Z){
    if(X<Dimension.at(0)&&Y<Dimension.at(1)&&Z<Dimension.at(2) && !(X<0) &&!(Y<0) &&!(Z<0)){
        return SimuNodes.at(X+Y*Dimension.at(0)+Z*Dimension.at(0)*Dimension.at(1));
    }
    else throw std::invalid_argument("Index Not in Mesh");
}

void PhaseSimulationMesh::showMeshProp(){
    std::cout<<"PhaseSimulationMesh Properties:"<<std::endl;
    std::cout<<"Mesh Size:\t\t"<<BoxX<<"\u0078"<<BoxY<<"\u0078"<<BoxY<<std::endl;
    std::cout<<"Number of Nodes:\t"<<getNumNodes()<<std::endl;
}

void PhaseSimulationMesh::showNodesProp(unsigned index){
    for(long i = 0; i < getNumNodes(); i++){
        // std::cout<<SimuNodes.at(i).getProperties(index)<<" ";
        if(i%BoxX==BoxX-1)std::cout<<std::endl;
        if(i%(BoxX*BoxY) == BoxX*BoxY-1)std::cout<<std::endl;
    }
}
#endif