#pragma once
#ifndef PHASE_SIMULATION_MESH
#define PHASE_SIMULATION_MESH

#include <iostream>
#include <vector>
#include "PhaseNode.hh"

enum WHICHDIM{X,Y,Z};
enum BOUNDCOND{PERIODIC,CONST,ADIABATIC};

class PhaseSimulationMesh{
    private:
        std::vector<int> Dimension{64, 64, 1};
        int& BoxX = Dimension.at(0);
        int& BoxY = Dimension.at(1);
        int& BoxZ = Dimension.at(2);
        std::vector<PhaseNode> SimuNodes;
    public:

        PhaseSimulationMesh(){} // initial with default size, without SimuNodes.at(i)s

        PhaseSimulationMesh(std::vector<int> SizeInfo, PhaseNode Nodes){ // initial with size and SimuNodes.at(i)s
            Dimension = SizeInfo;
            fillNodes(Nodes);
            bindBoundary(BOUNDCOND::PERIODIC);
        }
        PhaseSimulationMesh(PhaseNode Nodes){ // initial with SimuNodes.at(i)s and default size
            fillNodes(Nodes);
        }
        ~PhaseSimulationMesh(){};

/*************************************************************/

        void fillNodes(PhaseNode Nodes){ // fill mesh with SimuNodes.at(i)s 
            for(int i = 0; i < Dimension.at(0)*Dimension.at(1)*Dimension.at(2); i++){
                SimuNodes.push_back(Nodes);
            }
        }
        void bindBoundary(int whichBOUNDCOND){
            for(int i = 0; i < getNumNodes(); i++ ){
                SimuNodes.at(i).Forward = (i-BoxY < 0? &(SimuNodes.at(i)): &(SimuNodes.at(i-BoxY)) );
                SimuNodes.at(i).Backward = (i+BoxY >= getNumNodes()? &(SimuNodes.at(i)) : &(SimuNodes.at(i+BoxY)));
                SimuNodes.at(i).Left =(i-1 < 0? &(SimuNodes.at(i)) : &(SimuNodes.at(i-1)));
                SimuNodes.at(i).Right = (i+1 >= getNumNodes()? &(SimuNodes.at(i)) : &(SimuNodes.at(i+1)));
                SimuNodes.at(i).Down= (i-BoxX*BoxY < 0? &(SimuNodes.at(i)): &(SimuNodes.at(i-BoxX*BoxY)));
                SimuNodes.at(i).Up= (i+BoxX*BoxY >= getNumNodes()? &(SimuNodes.at(i)): &(SimuNodes.at(i+BoxX*BoxY)));

                if(whichBOUNDCOND == BOUNDCOND::PERIODIC){
                    if(i%BoxX == 0){
                        SimuNodes.at(i).Left = &(SimuNodes.at(i+(BoxX-1)));
                    }
                    if(i%BoxX == (BoxX-1)){
                        SimuNodes.at(i).Right = &(SimuNodes.at(i-(BoxX-1)));
                    }
                    if((i%(BoxX*BoxY))>=0 && (i%(BoxX*BoxY))<BoxX){
                        SimuNodes.at(i).Forward = &(SimuNodes.at(i+BoxX*(BoxY-1)));
                    }
                    if((i%(BoxX*BoxY))>=BoxX*(BoxY-1) && (i%(BoxX*BoxY))<BoxX*BoxY){
                        SimuNodes.at(i).Backward = &(SimuNodes.at(i-BoxX*(BoxY-1)));
                    }
                    if(i>=0 && i<BoxX*BoxY){
                        SimuNodes.at(i).Down = &(SimuNodes.at(i+BoxX*BoxY*(BoxZ-1)));
                    }
                    if(i>=BoxX*BoxY*(BoxZ-1) && i<BoxX*BoxY*BoxZ){
                        SimuNodes.at(i).Up = &(SimuNodes.at(i-BoxX*BoxY*(BoxZ-1)));
                    }
                }
            }
        }

        int getDim(const int which){
            return Dimension.at(which);
        }
        int getNumNodes(){ // return the number of SimuNodes.at(i)s in mesh
            return SimuNodes.size();
        }
        
        PhaseNode& findNode(int X, int Y, int Z){ // find SimuNodes.at(i) in the mesh according to the coordinates
            if(X<Dimension.at(0)&&Y<Dimension.at(1)&&Z<Dimension.at(2) && !(X<0) &&!(Y<0) &&!(Z<0)){
                return SimuNodes.at(X+Y*Dimension.at(0)+Z*Dimension.at(0)*Dimension.at(1));
            }
            else throw std::invalid_argument("Index Not in Mesh");
        }
        PhaseNode& findNode(int AbsCoord){
            if(AbsCoord<getNumNodes()){
                return SimuNodes.at(AbsCoord);
            }
            else throw std::invalid_argument("Index Not in Mesh");
        }

        void showMeshProp(); // show the basic information of the mesh
        void showNodesProp(unsigned which, unsigned index); // show one of the properties of the SimuNodes.at(i)s inside the mesh
};

void PhaseSimulationMesh::showMeshProp(){
    std::cout<<"PhaseSimulationMesh Properties:"<<std::endl;
    std::cout<<"Mesh Size:\t\t"<<BoxX<<"\u0078"<<BoxY<<"\u0078"<<BoxZ<<std::endl;
    std::cout<<"Number of Nodes:\t"<<getNumNodes()<<std::endl;
}

void PhaseSimulationMesh::showNodesProp(unsigned which, unsigned index){ //which para, index of para
    if(which == 0)std::cout<<"Concentration of "<<index<<" element\n";
    if(which == 1)std::cout<<"Order Parameter of "<<index<<" grain\n";
    for(long i = 0; i < getNumNodes(); i++){
        std::cout<<std::fixed<<std::setprecision(4)<<SimuNodes.at(i).getProperties(which).at(index)<<" ";
        if(i%BoxX==BoxX-1)std::cout<<"\n";
        if(i%(BoxX*BoxY) == BoxX*BoxY-1)std::cout<<"\n";
    }
    std::cout<<std::endl;
}
#endif