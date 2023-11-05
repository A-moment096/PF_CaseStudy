#pragma once
#ifndef PHASE_NODE
#define PHASE_NODE

#include <iostream>
#include <vector>
#include <random>
/***********************************************
This is PhaseNode Class, Every method is focused
only on one node (or say, point) in simulation
region (or, box).
************************************************/
std::mt19937 generator;


class PhaseNode{
    private:
        std::vector<long double> NodeProperties{298.15,0.0,0.0};
        long double& Temperature = NodeProperties.at(0);
        long double& Concentration = NodeProperties.at(1);
        long double& OrderParameter = NodeProperties.at(2);
        // long double PhaseFraction;
    public:
// Construct & Deconstruct Functions
        PhaseNode(){}; // accept default parameters
        PhaseNode(std::vector<long double> NodeInfo); //node from a vector
        PhaseNode(const PhaseNode &NewNode); // node from other node
        ~PhaseNode(){};
// Manipulate Methods
        long double getProperties(unsigned index); // return node properties 
        std::vector<long double> getProperties();

        void showNode();
        void updateNode(std::vector<long double> NewNode );

        void ConInitial_AveDis(long double Ave, long double Var); // parameter average distribution initialization
        // void ParaNormDisInitial(long double sigma, long double mu);
};


PhaseNode::PhaseNode(std::vector<long double> NodeInfo){
    Temperature = NodeInfo.at(0);
    Concentration = NodeInfo.at(1);
    OrderParameter = NodeInfo.at(1);
}

PhaseNode::PhaseNode( const PhaseNode &NewNode){
    Temperature = NewNode.Temperature;
    Concentration = NewNode.Concentration;
    OrderParameter = NewNode.OrderParameter;
}

void PhaseNode::ConInitial_AveDis(long double Ave, long double Var){
    Concentration = Ave +Var-2*double(rand()%(1000*int(Var)))/1000;
}

void PhaseNode::updateNode(std::vector<long double> NewNode){
    NodeProperties = NewNode;
}

void PhaseNode::showNode(){
    std::cout<<"Node Information:"<<std::endl;
    std::cout<<"Temperature:\t\t"<<Temperature<<std::endl;
    std::cout<<"Concentration:\t\t"<<Concentration<<std::endl;
    std::cout<<"OrderParameter:\t\t"<<OrderParameter<<std::endl;
}

std::vector<long double> PhaseNode::getProperties(){
    return NodeProperties;
}

long double PhaseNode::getProperties(unsigned index){
    return NodeProperties.at(index);
}

//void PhaseNode::ParaNormDisInitial(long double mu, long double sigma){
//    while(true){
//        std::normal_distribution<long double> randnum(mu,sigma);
//        long double r = randnum(generator);
//        if(r > mu-3*sigma && r < mu+3*sigma){
//            Concentration = r;
//            break;
//        }
//    }
//}


#endif