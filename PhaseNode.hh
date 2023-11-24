#pragma once
#ifndef PHASE_NODE
#define PHASE_NODE

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include "PhaseEntry.hh"

class PhaseNode{
    public:
        unsigned Num_Ent;
        std::vector<PhaseEntry> Entrys;

/*************************************************************/

// Construct & Deconstruct Functions
        PhaseNode(){
            Entrys.push_back(Def_PhsEnt);
            Num_Ent = 1;
            (*this)(0).Index = 0;
        } //Accept Default Parameters
        
        PhaseNode(std::vector<PhaseEntry> EntryList){
            Entrys = EntryList;
            Num_Ent = EntryList.size();
            updateIndex();
        }

        PhaseNode(const PhaseNode &_node){
            Num_Ent = _node.Num_Ent;
            Entrys = _node.Entrys;
        }

        PhaseNode& operator= (const PhaseNode &_node){
            Num_Ent = _node.Num_Ent;
            Entrys = _node.Entrys;
            return *this;
        }

        PhaseEntry& operator() (const unsigned _Index){
            if(_Index<Num_Ent)
            return Entrys.at(_Index);
            else throw std::out_of_range("No such entry");
        }

        ~PhaseNode(){
            Entrys.clear();
        }

        /*************************************************************/

        std::vector<double> getVal(){
            std::vector<double> result;
            for(auto ent : Entrys){
                result.push_back(ent.PhaseFrac);
            }
            return result;
        }

        std::vector<double> getLap(){
            std::vector<double> result(Num_Ent);
            for(auto ent : Entrys){
                result.push_back(ent.Lap);
            }
            return result;
        }

        std::vector<double> getGrad(){
            std::vector<double> result(Num_Ent);
            for(auto ent : Entrys){
                result.push_back(ent.Grad);
            }
            return result;
        }
        double getVal(unsigned _Index){
            if(_Index < Num_Ent)
            return (*this)(_Index).PhaseFrac;
            else throw std::out_of_range("No Such Index");
        }

        double getLap(unsigned _Index){
            if(_Index < Num_Ent)
            return (*this)(_Index).Lap;
            else throw std::out_of_range("No Such Index");
        }

        double getGrad(unsigned _Index){
            if(_Index < Num_Ent)
            return (*this)(_Index).Grad;
            else throw std::out_of_range("No Such Index");
        }

        /*****************************************************************/

        void addEntry(int num){
            for (int i = 0; i < num; ++i)
            {
                Entrys.push_back(Def_PhsEnt);
            }
            Num_Ent = Entrys.size();
            updateIndex();
        }

        void deletEntry(unsigned _Index){
            if(Num_Ent <= 1) throw std::out_of_range("Last entry");
            else
                for(auto &ent : Entrys){
                    if(ent.Index == _Index){
                        std::vector<PhaseEntry>::iterator delPos = std::find(Entrys.begin(),Entrys.end(),ent);
                        Entrys.erase(delPos);
                        return;
                    }
                }
            throw std::out_of_range("Not in entry list");
        }

        void updateVal(unsigned _Index, double _PhsFrac){
            for(auto &ent : Entrys){
                if(ent.Index == _Index){
                    ent.PhaseFrac = _PhsFrac;
                    return;
                }
            }
            throw std::out_of_range("Not in entry list");
        }

        void updateLap(unsigned _Index, double _Lap){
            for(auto &ent : Entrys){
                if(ent.Index == _Index){
                    ent.Lap = _Lap;
                    return;
                }
            }
            throw std::out_of_range("Not in entry list");
        }

        void updateGrad(unsigned _Index, double _Grad){
            for(auto &ent : Entrys){
                if(ent.Index == _Index){
                    ent.Grad = _Grad;
                    return;
                }
            }
            throw std::out_of_range("Not in entry list");
        }

        /*****************************************************************/

        void updateIndex(){
            for(unsigned i = 0; i < Num_Ent; i++){
                (*this)(i).Index = (i);
            }
        }

        double sumPhsFrac(){
            double result = 0;
            for(auto ent : Entrys)result += ent.PhaseFrac;
            return result;  
        }

        double sumPhsFrac2(){
            double result = 0;
            for(auto ent : Entrys)result += ent.PhaseFrac*ent.PhaseFrac;
            return result;
        }

        double sumPhsFrac3(){
            double result = 0;
            for(auto ent : Entrys)result += ent.PhaseFrac*ent.PhaseFrac*ent.PhaseFrac;
            return result;
        }


}Def_PhsNode;

#endif
