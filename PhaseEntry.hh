#pragma once
#ifndef PHASE_ENTRY
#define PHASE_ENTRY

#include <iostream>
#include <vector>

class PhaseEntry{
    public:
        double PhaseFrac = 0;
        double Lap = 0;
        double Grad = 0;
        unsigned Index;
        
/*************************************************************/

        PhaseEntry(){
            Index = 0;
        }

        PhaseEntry(const PhaseEntry& _Entry){
            PhaseFrac = _Entry.PhaseFrac;
            Index = _Entry.Index;
            Lap = _Entry.Lap;
            Grad = _Entry.Grad;
        }

        PhaseEntry& operator= (const PhaseEntry& _Entry){
            PhaseFrac = _Entry.PhaseFrac;
            Index = _Entry.Index;
            Lap = _Entry.Lap;
            Grad = _Entry.Grad;
            return *this;
        }

        ~PhaseEntry(){
        }

        
/**/    bool operator== (const PhaseEntry& _Entry){
            bool result;
            ((PhaseFrac == _Entry.PhaseFrac)&&(Index == _Entry.Index)) ? (result = true) : (result = false);
            return result;
        }

/*************************************************************/


}Def_PhsEnt;

#endif