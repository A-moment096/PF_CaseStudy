#pragma once
#ifndef PHASE_NODE
#define PHASE_NODE

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include "PhaseEntry.hh"

class PhaseNode{
    private:
        unsigned Num_PhsFrc;
        std::vector<PhaseEntry> Entrys;
    public:
// Construct & Deconstruct Functions
        PhaseNode(){
            Entrys.push_back(Def_PhsEnt);
            Num_PhsFrc = 1;
        } //Accept Default Parameters
        
        PhaseNode(std::vector<PhaseEntry> EntryList){
            Entrys = EntryList;
            Num_PhsFrc = EntryList.size();
            for(unsigned i = 0; i < Num_PhsFrc; i++){
                Entrys.at(i).setindex(i);
            }
        }

        PhaseNode(const PhaseNode &_node){
            Num_PhsFrc = _node.Num_PhsFrc;
            Entrys = _node.Entrys;
        }

        PhaseNode& operator= (const PhaseNode &_node){
            Num_PhsFrc = _node.Num_PhsFrc;
            Entrys = _node.Entrys;
            return *this;
        }

        ~PhaseNode(){
            Entrys.clear();
        }

        unsigned getNums(){
            return Num_PhsFrc;
        }

        std::vector<double> getPhsFrac(){
            std::vector<double> result;
            for(auto ent : Entrys){
                result.push_back(ent.getPhsFrac());
            }
            return result;
        }

        double getPhsFrac(unsigned index){
            if( index < Num_PhsFrc){
                return Entrys.at(index).getPhsFrac();
            }
            else throw std::out_of_range("Not in entry list");
        }

        std::vector<unsigned> getindex(){
            std::vector<unsigned> result;
            for(auto ent : Entrys){
                result.push_back(ent.getindex());
            }
            return result;
        }

        PhaseEntry getEntry(unsigned index){
            for(auto ent : Entrys){
                if(ent.getindex() == index){
                    return ent;
                }
            }
            throw std::out_of_range("Not in entry list");
        }


        void addEntry(PhaseEntry _Entry){
            Entrys.push_back(_Entry);
        }

        void deletEntry(unsigned index){
            if(Num_PhsFrc <= 1) throw std::out_of_range("Last entry");
            else
                for(auto &ent : Entrys){
                    if(ent.getindex() == index){
                        std::vector<PhaseEntry>::iterator delPos = std::find(Entrys.begin(),Entrys.end(),ent);
                        Entrys.erase(delPos);
                        return;
                    }
                }
            throw std::out_of_range("Not in entry list");
        }

        void updateEntry(unsigned index, double _PhsFrac){
            for(auto &ent : Entrys){
                if(ent.getindex() == index){
                    ent.setPhsFrac(_PhsFrac);
                    return;
                }
            }
            throw std::out_of_range("Not in entry list");
        }
        
        void updateEntry(double _phs){
            if(Num_PhsFrc == 1){
                Entrys.at(0).setPhsFrac(_phs);
                return;
            }
            else throw std::invalid_argument("Not only one element");
        }

}Def_PhsNode;

#endif
