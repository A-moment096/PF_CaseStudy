#pragma once
#ifndef CON_NODE
#define CON_NODE

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include "ConEntry.hh"

class ConNode{
    public:
        unsigned Num_Ent;
        std::vector<ConEntry> Entrys;

/*************************************************************/

// Construct & Deconstruct Functions
        ConNode(){
            Entrys.push_back(Def_ConEnt);
            Num_Ent = 1;
            (*this)(0).Index  = 0;
        } //Accept Default Parameters
        
        ConNode(std::vector<ConEntry> EntryList){
            Entrys = EntryList;
            Num_Ent = EntryList.size();
            updateIndex();
        }

        ConNode(std::vector<ELEMENT> ElementList){
            Num_Ent = ElementList.size();
            for (unsigned i = 0; i < Num_Ent; i++){
                Entrys.push_back(Def_ConEnt);
                (*this)(i).Element = (ElementList.at(i));
                (*this)(i).Index  = i;
            }
        }

        ConNode(const ConNode& _node){
            Num_Ent = _node.Num_Ent;
            Entrys = _node.Entrys;
        }

        ConNode& operator= (const ConNode& _node){
            Num_Ent = _node.Num_Ent;
            Entrys = _node.Entrys;
            return *this;
        }
                
        ConEntry& operator() (const unsigned _Index){
            if(_Index<Num_Ent)
            return Entrys.at(_Index);
            else throw std::out_of_range("No such entry");
        }

        ~ConNode(){
            Entrys.clear();
        }

        /*************************************************************/

        std::vector<double> getVal(){
            std::vector<double> result;
            for(auto ent : Entrys){
                result.push_back(ent.Concentration);
            }
            return result;
        }

        std::vector<double> getLap(){
            std::vector<double> result;
            for(auto ent : Entrys){
                result.push_back(ent.Lap);
            }
            return result;
        }

        std::vector<double> getGrad(){
            std::vector<double> result;
            for(auto ent : Entrys){
                result.push_back(ent.Grad);
            }
            return result;
        }

        double getVal(unsigned _Index){
            if(_Index < Num_Ent)
            return (*this)(_Index).Concentration;
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

        /*************************************************************/

        void addEntry(unsigned num){
            for (unsigned i = 0; i < num; ++i)
            {
                Entrys.push_back(Def_ConEnt);
            }
            Num_Ent = Entrys.size();
            updateIndex();
        }

        void deletEntry(ELEMENT delElemnt){
            if(Num_Ent <= 1) throw std::out_of_range("Last entry");
            else
                for(auto &ent : Entrys){
                    if(ent.Element == delElemnt){
                        std::vector<ConEntry>::iterator delPos = std::find(Entrys.begin(),Entrys.end(),ent);
                        Entrys.erase(delPos);
                        return;
                    }
                }
            throw std::invalid_argument("No such element");
        }

        void deletEntry(unsigned _Index){
            if(Num_Ent <= 1) throw std::out_of_range("Last entry");
            else
                for(auto &ent : Entrys){
                    if(ent.Index == _Index){
                        std::vector<ConEntry>::iterator delPos = std::find(Entrys.begin(),Entrys.end(),ent);
                        Entrys.erase(delPos);
                        return;
                    }
                }
            throw std::out_of_range("Not in entry list");
        }

        void updateVal(ELEMENT _elemnt,double _con){
            for(auto &ent : Entrys){
                if(ent.Element == _elemnt){
                    ent.Concentration = _con;
                    return;
                }
            }
            throw std::invalid_argument("No such element");
        }

        void updateVal(unsigned _Index,double _con){
            for(auto &ent : Entrys){
                if(ent.Index == _Index){
                    ent.Concentration = _con;
                    return;
                }
            }
            throw std::invalid_argument("No such element");
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
    
        /*************************************************************/

        void updateIndex(){
            for(unsigned i = 0; i < Num_Ent; i++){
                (*this)(i).Index = i;
            }
        }
}Def_ConNode;

#endif