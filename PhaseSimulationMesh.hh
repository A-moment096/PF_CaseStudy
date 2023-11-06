#pragma once
#ifndef PHASE_SIMULATION_MESH
#define PHASE_SIMULATION_MESH

#include <iostream>
#include <vector>

class PhaseSimulationMesh{
    private:
        std::vector<long int> Dimension{64, 64, 1};
        long int& BoxX = Dimension.at(0);
        long int& BoxY = Dimension.at(1);
        long int& BoxZ = Dimension.at(2);

        // std::vector<PhaseNode> SimuNodes;

    public:
        PhaseSimulationMesh(){}; // initial with default size, without nodes
        PhaseSimulationMesh(const PhaseSimulationMesh &NewMesh );
        PhaseSimulationMesh(std::vector<long int> MeshInfo);
        ~PhaseSimulationMesh(){};
        PhaseSimulationMesh& operator= (const PhaseSimulationMesh &NewMesh);

/**/    void showMeshProp(); // show the basic information of the mesh

        // long getNumNodes(); // return the number of nodes in mesh
        // void fillNodes(PhaseNode Nodes); // fill mesh with nodes 

/**/    // void Ave_Dis_Noise(long double Ave, long double Var); // add noise to the nodes (concentration) with average distribution
        // PhaseNode& findNode(long X, long Y, long Z); // find node in the mesh according to the coordinates

        // void showNodesProp(unsigned index); // show one of the properties of the nodes inside the mesh
} DefaultMesh;

PhaseSimulationMesh::PhaseSimulationMesh(const PhaseSimulationMesh &NewMesh ){
    Dimension = NewMesh.Dimension;
}

PhaseSimulationMesh::PhaseSimulationMesh(std::vector<long> MeshInfo){
    Dimension = MeshInfo;
}

PhaseSimulationMesh& PhaseSimulationMesh::operator= (const PhaseSimulationMesh &NewMesh){
    Dimension = NewMesh.Dimension;
    return *this;
}


void PhaseSimulationMesh::showMeshProp(){
    std::cout<<"PhaseSimulationMesh Properties:"<<std::endl;
    std::cout<<"Mesh Size:\t\t"<<BoxX<<"\u0078"<<BoxY<<"\u0078"<<BoxZ<<std::endl;
    std::cout<<"Number of Nodes:\t"<<BoxX*BoxY*BoxZ<<std::endl;
}


// long PhaseSimulationMesh::getNumNodes(){
//     return SimuNodes.size();
// }
// 
// void PhaseSimulationMesh::fillNodes(PhaseNode Nodes){
//     for(long i = 0; i < BoxX*BoxY*BoxZ; i++){
//         SimuNodes.push_back(Nodes);
//     }
// }
// 
// void PhaseSimulationMesh::Ave_Dis_Noise(long double Ave, long double Var){
//     for(auto &Nodes : SimuNodes){
//         Nodes.ConInitial_AveDis(Ave,Var);
//     }
// }
// 
// 
// PhaseNode& PhaseSimulationMesh::findNode(long X, long Y, long Z){
//     if(X<BoxX&&Y<BoxY&&Z<BoxY && !(X<0) &&!(Y<0) &&!(Z<0)){
//         return SimuNodes.at(X+Y*BoxX+Z*BoxX*BoxY);
//     }
//     else throw std::invalid_argument("Index Not in Mesh");
// }
// void PhaseSimulationMesh::showNodesProp(unsigned index){
//     for(long i = 0; i < getNumNodes(); i++){
//         std::cout<<SimuNodes.at(i).getProperties(index)<<" ";
//         if(i%BoxX==BoxX-1)std::cout<<std::endl;
//         if(i%(BoxX*BoxY) == BoxX*BoxY-1)std::cout<<std::endl;
//     }
// }

#endif