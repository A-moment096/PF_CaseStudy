#pragma once
#ifndef PHASE_NODE
#define PHASE_NODE

#include <iostream>
#include <iomanip>
#include <vector>

#include "PhaseEntry.hh"
#include "ConEntry.hh"
/***********************************************
This is PhaseNode Class, Every method is focused
only on one node (or say, point) in simulation
region (or, box).
************************************************/

enum WHICHPARA {CON,GRAIN,CUSTOM=99};

class PhaseNode{
    private:
        unsigned Num_Elemt = 1;
        unsigned Num_Grain = 1;
        std::vector<ConEntry> ConVect;
        std::vector<PhaseEntry> PhsVect;
        double CustomValue=0;
        double Temperature = 298.15;
    public:
        PhaseNode* Up = this;
        PhaseNode* Down = this;
        PhaseNode* Forward = this;
        PhaseNode* Backward = this;
        PhaseNode* Left = this;
        PhaseNode* Right = this;
// Construct & Deconstruct Functions
        PhaseNode(){
            PhsVect.push_back(DefaultPhsEnt);
            ConVect.push_back(DefaultConEnt);
        }; //Accept Default Parameters
    
        PhaseNode(unsigned n_elemt, unsigned n_grain):Num_Elemt(n_elemt),Num_Grain(n_grain){
            for (int i = 0; i < Num_Elemt; i++)
            {
                ConVect.push_back(DefaultConEnt);
                ConVect.at(i).setindex(i);
            }
            for (int j = 0; j < Num_Grain; j++)
            {
                PhsVect.push_back(DefaultPhsEnt);
                PhsVect.at(j).setindex(j);
            }
        }

        PhaseNode(unsigned n_grain):PhaseNode(1,n_grain){}

        PhaseNode(const PhaseNode &NewNode){
            ConVect = NewNode.ConVect;
            PhsVect = NewNode.PhsVect;
            Temperature = NewNode.Temperature;
        } //Node from Other Node
        PhaseNode& operator= (const PhaseNode& NewNode){
            ConVect = NewNode.ConVect;
            PhsVect = NewNode.PhsVect;
            Temperature = NewNode.Temperature;
            return *this;
        }
        ~PhaseNode(){};
// Manipulate Methods
        
        unsigned getNums(int which){
            if(which == WHICHPARA::CON){
                return Num_Elemt;
            }
            if(which == WHICHPARA::GRAIN){
                return Num_Grain;
            }
        }

        std::vector<double> getProperties(int which){
            
            if(which == WHICHPARA::CON){ // concentration
            std::vector<double> temp(Num_Elemt,0.0);
                for (int i = 0; i < Num_Elemt; i++)
                {
                    temp.at(i) = (ConVect.at(i).getCon());
                }
                return temp;
            }

            if(which == WHICHPARA::GRAIN){ // grain order parameter
            std::vector<double> temp(Num_Grain,0.0);
                for (int i = 0; i < Num_Grain; i++)
                {
                    temp.at(i) = (PhsVect.at(i).getOrderPara());
                }
                return temp;
            }

            if(which == WHICHPARA::CUSTOM){
            std::vector<double> temp(1);
                temp.at(0) = (CustomValue);
                return temp;
            }
        }

        void updateNode(unsigned which, unsigned index, double value){
            if(which == WHICHPARA::CON){
                ConVect.at(index).setCon(value);
            }
            if(which == WHICHPARA::GRAIN){
                PhsVect.at(index).setOrderPara(value);
            }
            if(which == WHICHPARA::CUSTOM){
                CustomValue = value;
            }
        }

        void showNode();
        // void ConInitial_AveDis(double Ave, double Var); //Parameter Average Distribution Initialization
};

/*************************************************************/
void PhaseNode::showNode(){
    std::cout<<"Node Information:\n";
    std::cout<<"Temperature:\t\t"<<Temperature<<"\n";
    std::cout<<"OrdParaIndex:\tOrderParameter:\t\tConIndex:\tConcentration:\n";
    for(int j = 0; j < Num_Grain; j++){
        for(int i = 0; i < Num_Elemt ; i++){
            std::cout<<PhsVect.at(j).getindex()<<"\t\t"<<std::fixed<<std::setprecision(6)<<PhsVect.at(j).getOrderPara()<<"\t\t";
            std::cout<<ConVect.at(i).getindex()<<"\t\t"<<std::fixed<<std::setprecision(6)<<ConVect.at(i).getCon()<<"\n";
        }
        std::cout<<"\n";
    }
    std::cout<<std::endl;
}


#endif