#pragma once
#ifndef PHASE_ENTRY
#define PHASE_ENTRY

#include <iostream>
#include <vector>

class PhaseEntry{
    private:
    double OrderParameter = 0;
    unsigned index;
    public:
    
    PhaseEntry(){};
    ~PhaseEntry(){};
    PhaseEntry& operator= (const PhaseEntry& newEntry){
            OrderParameter = newEntry.OrderParameter;
            index = newEntry.index;
            return *this;
        }

    double getOrderPara(){
        return OrderParameter;
    }
    int getindex(){
        return index;
    }
    void setOrderPara(const double newOrderPara ){
        OrderParameter = newOrderPara;
    }
    void setindex(const int _index){
        index = _index;
    }


}DefaultPhsEnt;

#endif