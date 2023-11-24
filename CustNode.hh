#pragma once

#ifndef CUST_NODE
#define CUST_NODE

#include <iostream>
#include <iomanip>
#include <vector>

class CustNode
{
public:
    std::vector<double> CustVal;
    std::vector<double> CustLap;
    std::vector<double> CustGrad;
    int Num;

    CustNode():CustNode(1){}

    CustNode(int _num){
        Num = _num;
        CustVal.reserve(Num);
        CustLap.reserve(Num);
        CustGrad.reserve(Num);
    }

    ~CustNode(){
        CustVal.clear();
        CustLap.clear();
        CustGrad.clear();
    }

    void updateLap(int _index, double _val){
        CustLap.at(_index) = _val;
    }

    void addEntry(int num){
        Num += num;
        CustVal.reserve(Num);
        CustLap.reserve(Num);
        CustGrad.reserve(Num);
    }

    std::vector<double> getVal(){
        return CustVal;
    }

    double getVal(int _index){
        return CustVal.at(_index);
    }

    double getLap(int _index){
        return CustLap.at(_index);
    }

    double getGrad(int _index){
        return CustGrad.at(_index);
    }

}Def_CustNode;


#endif