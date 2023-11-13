#pragma once
#ifndef CON_ENTRY
#define CON_ENTRY

#include <iostream>
#include <vector>

class ConEntry{
    private:
    double Concentration = 0;
    unsigned index;
    public:
    
    ConEntry(){};
    ~ConEntry(){};
    ConEntry& operator= (const ConEntry& newEntry){
            Concentration = newEntry.Concentration;
            return *this;
        }

    double getCon(){
        return Concentration;
    }
    int getindex(){
        return index;
    }
    void setCon(const double newOrderPara ){
        Concentration = newOrderPara;
    }
    void setindex(const int _index){
        index = _index;
    }

}DefaultConEnt;

#endif