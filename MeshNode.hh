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

enum WHICHPARA {CON,PHSFRAC,CUSTOM=99};
// enum WHICHDIR {UP,DOWN,FORWARD,BACKWARD,LEFT,RIGHT};

class MeshNode{
    public:
        
        double Temperature = 298.15;

        PhaseNode Phs_Node;
        ConNode Con_Node;
        
        std::vector<double> Custom_Value;
        double CustLap = 0.0;

        MeshNode* Up = nullptr;
        MeshNode* Down = nullptr;
        MeshNode* Forward = nullptr;
        MeshNode* Backward = nullptr;
        MeshNode* Left = nullptr;
        MeshNode* Right = nullptr;

/*************************************************************/

     // Construct & Deconstruct Functions
        MeshNode(){
            Phs_Node = Def_PhsNode;
            Con_Node = Def_ConNode;
            Custom_Value.reserve(10);
        }; //Accept Default Parameters
    
        MeshNode(PhaseNode _phs_node, ConNode _con_node){
            Phs_Node = _phs_node;
            Con_Node = _con_node;
            Custom_Value.reserve(10);
        }

        MeshNode(PhaseNode _phs_node):MeshNode(_phs_node, Def_ConNode){}
        
        MeshNode(ConNode _con_node):MeshNode(Def_PhsNode,_con_node){}
        
        ~MeshNode(){
        Up = nullptr;
        Down = nullptr;
        Forward = nullptr;
        Backward = nullptr;
        Left = nullptr;
        Right = nullptr;
        // Custom_Value.clear();
        };

        /*************************************************************/
        // Manipulate Methods

        unsigned getNum_Ent(WHICHPARA whichpara){
            switch (whichpara)
            {
            case WHICHPARA::CON :
                return Con_Node.Num_Ent;
                break;
            case WHICHPARA::PHSFRAC :
                return Phs_Node.Num_Ent;
                break;
            case WHICHPARA::CUSTOM :
                return 1;
                break;
            default:
                throw std::invalid_argument("No such element");
                return 0;
                break;
            }
        }

        std::vector<double> getProp(WHICHPARA whichpara){
            switch (whichpara)
            {
            case WHICHPARA::CON :
                return Con_Node.getCon();
                break;
            case WHICHPARA::PHSFRAC :
                return Phs_Node.getPhsFrac();
                break;
            case WHICHPARA::CUSTOM :
                return {Custom_Value};
                break;
            default:
                break;
            }
            return {};
        }

        std::vector<double> getLap(unsigned whichpara){
            switch (whichpara)
            {
            case WHICHPARA::CON :
                return Con_Node.getLap();
                break;
            case WHICHPARA::PHSFRAC :
                return Phs_Node.getLap();
                break;
            default:
                break;
            }
            return {};
        }

         std::vector<double> getGrad(unsigned whichpara){
            switch (whichpara)
            {
            case WHICHPARA::CON :
                return Con_Node.getGrad();
                break;
            case WHICHPARA::PHSFRAC :
                return Phs_Node.getGrad();
                break;
            default:
                break;
            }
            return {};
        }
        
        /*************************************************************/

        double sumPhsFrac(){
            return Phs_Node.sumPhsFrac();
        }

        double sumPhsFrac2(){
            return Phs_Node.sumPhsFrac2();
        }
        
        double sumPhsFrac3(){
            return Phs_Node.sumPhsFrac3();
        }
        

        /*************************************************************/

        void showNode();
}Def_Node;

/*************************************************************/
void MeshNode::showNode(){
    std::cout<<"Node Information:\n";
    std::cout<<"Temperature:\t\t"<<Temperature<<"\n";
    std::cout<<"Phase Index:\tPhase Fraction:\t\tElement:\tConcentration:\n";
    for(int i = 0; i < Phs_Node.Num_Ent; i++){
        for(int j = 0; j < Con_Node.Num_Ent; j++){
            std::cout<<Phs_Node.Entrys.at(i).index<<"\t\t"<<std::fixed<<std::setprecision(6)<<Phs_Node.getPhsFrac().at(i)<<"\t\t";
            std::cout<<Con_Node.Entrys.at(j).Element<<"\t\t"<<std::fixed<<std::setprecision(6)<<Con_Node.getCon().at(j)<<"\n";
        }
        std::cout<<"\n";
    }
    std::cout<<std::endl;
}


#endif

