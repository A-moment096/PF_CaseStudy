#pragma once

#ifndef CUST_NODE_HH
#define CUST_NODE_HH

#include <iostream>
#include <iomanip>
#include <vector>

class CustNode{
    public:
    std::vector<double> CustVal;
    std::vector<double> CustLap;
    std::vector<double> CustGrad;
    int Num;

    CustNode():CustNode(1){}

    CustNode(int _num){
        addEntry(_num);
    }

    ~CustNode(){
        CustVal.clear();
        CustLap.clear();
        CustGrad.clear();
    }

    void updateLap(int _Index, double _val){
        CustLap.at(_Index) = _val;
    }

    void addEntry(int _num){
        Num += _num;
        for (int i = 0; i<_num; ++i){
            CustVal.push_back(0.0);
            CustLap.push_back(0.0);
            CustGrad.push_back(0.0);
        }
    }

    std::vector<double> getVal(){
        return CustVal;
    }

    double getVal(int _Index){
        return CustVal.at(_Index);
    }

    double getLap(int _Index){
        return CustLap.at(_Index);
    }

    double getGrad(int _Index){
        return CustGrad.at(_Index);
    }

}Def_CustNode;


#endif