#pragma once
#ifndef PHASE_ENTRY
#define PHASE_ENTRY

#include <iostream>
#include <vector>

class PhaseEntry{
    private:
        double PhaseFrac = 0;
        unsigned index;
    public:

        PhaseEntry(){
            index = 0;
        }

        PhaseEntry(const PhaseEntry& _Entry){
            PhaseFrac = _Entry.PhaseFrac;
            index = _Entry.index;
        }

        PhaseEntry& operator= (const PhaseEntry& _Entry){
            PhaseFrac = _Entry.PhaseFrac;
            index = _Entry.index;
            return *this;
        }

        ~PhaseEntry(){
            std::cout<<"phsent deleted\n";
        }

/***    **********************************************************/
        
        bool operator== (const PhaseEntry& _Entry){
            bool result;
            ((PhaseFrac == _Entry.PhaseFrac)&&(index == _Entry.index)) ? (result = true) : (result = false);
            return result;
        }

        double getPhsFrac(){
            return PhaseFrac;
        }
        int getindex(){
            return index;
        }
        void setPhsFrac(const double _PhsFrc ){
            PhaseFrac = _PhsFrc;
        }
        void setindex(const int _index){
            index = _index;
        }


}Def_PhsEnt;

#endif