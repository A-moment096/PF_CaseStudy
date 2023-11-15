#pragma once
#ifndef MESH_NODE
#define MESH_NODE

#include <iostream>
#include <iomanip>
#include <vector>

#include "ConNode.hh"
#include "PhaseNode.hh"
/***********************************************
This is MeshNode Class, Every method is focused
only on one node (or say, point) in simulation
region (or, box).
************************************************/


class MeshNode{
    private:
        double Custom_Value=0.0;
        double Temperature = 298.15;
        PhaseNode Phs_Node;
        ConNode Con_Node;
    public:
        MeshNode* Up = this;
        MeshNode* Down = this;
        MeshNode* Forward = this;
        MeshNode* Backward = this;
        MeshNode* Left = this;
        MeshNode* Right = this;
     // Construct & Deconstruct Functions
        MeshNode(){
            Phs_Node = Def_PhsNode;
            Con_Node = Def_ConNode;
        }; //Accept Default Parameters
    
        MeshNode(PhaseNode _phs_node, ConNode _con_node){
            Phs_Node = _phs_node;
            Con_Node = _con_node;
        }

        MeshNode(PhaseNode _phs_node):MeshNode(_phs_node, Def_ConNode){}
        
        MeshNode(ConNode _con_node):MeshNode(Def_PhsNode,_con_node){}
        
        ~MeshNode(){};
     // Manipulate Methods

        void showNode();
}Def_Node;

/*************************************************************/
void MeshNode::showNode(){
    std::cout<<"Node Information:\n";
    std::cout<<"Temperature:\t\t"<<Temperature<<"\n";
    std::cout<<"OrdParaIndex:\tOrderParameter:\t\tElement:\tConcentration:\n";
    for(int i = 0; i < Phs_Node.getNums(); i++){
        for(int j = 0; j < Con_Node.getNums(); j++){
            std::cout<<Phs_Node.getindex().at(i)<<"\t\t"<<std::fixed<<std::setprecision(6)<<Phs_Node.getPhsFrac(i)<<"\t\t";
            std::cout<<Con_Node.getElement().at(j)<<"\t\t"<<std::fixed<<std::setprecision(6)<<Con_Node.getCon().at(j)<<"\n";
        }
        std::cout<<"\n";
    }
    std::cout<<std::endl;
}


#endif