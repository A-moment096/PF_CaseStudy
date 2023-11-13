#pragma once
#ifndef PHASE_NODE
#define PHASE_NODE

#include <iostream>
#include <iomanip>
#include <vector>
#include <random>

#include "PhaseEntry.hh"
#include "ConEntry.hh"
/***********************************************
This is PhaseNode Class, Every method is focused
only on one node (or say, point) in simulation
region (or, box).
************************************************/

class PhaseNode{
    private:
        unsigned Num_Grain = 1;
        unsigned Num_Elemt = 1;
        std::vector<ConEntry> ConVect;
        std::vector<PhaseEntry> PhsVect;
        double Temperature = 298.15;
    public:
// Construct & Deconstruct Functions
        PhaseNode(){
            PhsVect.push_back(DefaultPhsEnt);
            ConVect.push_back(DefaultConEnt);
        }; //Accept Default Parameters
    
        PhaseNode(unsigned n_grain, unsigned n_elemt):Num_Elemt(n_elemt),Num_Grain(n_grain){
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

        PhaseNode(unsigned n_grain):PhaseNode(n_grain,1){}

        PhaseNode(std::vector<double> NodeInfo); //Node from a vector
        PhaseNode(const PhaseNode &NewNode); //Node from Other Node
        ~PhaseNode(){};
        PhaseNode& operator= (const PhaseNode& NewNode){
            ConVect = NewNode.ConVect;
            PhsVect = NewNode.PhsVect;
            Temperature = NewNode.Temperature;
            return *this;
        }
// Manipulate Methods
        void showNode();
        
        void updateNode(double Node){};

        // void ConInitial_AveDis(double Ave, double Var); //Parameter Average Distribution Initialization
        // void ParaNormDisInitial(double sigma, double mu);
};


// PhaseNode::PhaseNode(std::vector<double> NodeInfo){
//     Temperature = NodeInfo.at(0);
//     Concentration = NodeInfo.at(1);
//     OrderParameter = NodeInfo.at(1);
// }

inline PhaseNode::PhaseNode( const PhaseNode &NewNode){
    ConVect = NewNode.ConVect;
    PhsVect = NewNode.PhsVect;
    Temperature = NewNode.Temperature;
}
PhaseNode& PhaseNode::operator= (const PhaseNode NewNode){
    Temperature = NewNode.Temperature;
    Concentration = NewNode.Concentration;
    OrderParameter = NewNode.OrderParameter;
    return *this;
}


// void PhaseNode::ConInitial_AveDis(double Ave, double Var){
//     Concentration = Ave +Var-2*double(rand()%(1000*int(Var)))/1000;
// }

void PhaseNode::updateNode(unsigned which, long double NewProp){
    NodeProperties.at(which) = NewProp;
}
/*************************************************************/
void PhaseNode::showNode(){
    std::cout<<"Node Information:\n";
    std::cout<<"Temperature:\t\t"<<Temperature<<"\n";
    std::cout<<"OrdParaIndex:\tOrderParameter:\t\tConIndex:\tConcentration:\n";
    for(int j = 0; j < Num_Grain; j++){
        for(int i = 0; i < Num_Elemt ; i++){
            PhsVect.at(j).setOrderPara(double(j)/12+double(i)*0.3);
            ConVect.at(i).setCon(double(i+5)/0.453-1-double(j));
            std::cout<<PhsVect.at(j).getindex()<<"\t\t"<<std::fixed<<std::setprecision(6)<<PhsVect.at(j).getOrderPara()<<"\t\t";
            std::cout<<ConVect.at(i).getindex()<<"\t\t"<<std::fixed<<std::setprecision(6)<<ConVect.at(i).getCon()<<"\n";
        }
        std::cout<<"\n";
    }
    std::cout<<std::endl;
}

//void PhaseNode::ParaNormDisInitial(double mu, double sigma){
//    while(true){
//        std::normal_distribution<double> randnum(mu,sigma);
//        double r = randnum(generator);
//        if(r > mu-3*sigma && r < mu+3*sigma){
//            Concentration = r;
//            break;
//        }
//    }
//}


#endif