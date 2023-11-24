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
            Entrys.at(0).index  = 0;
        } //Accept Default Parameters
        
        ConNode(std::vector<ConEntry> EntryList){
            Entrys = EntryList;
            Num_Ent = EntryList.size();
            update_index();
        }

        ConNode(std::vector<ELEMENT> ElementList){
            Num_Ent = ElementList.size();
            for (unsigned i = 0; i < Num_Ent; i++){
                Entrys.push_back(Def_ConEnt);
                Entrys.at(i).Element = (ElementList.at(i));
                Entrys.at(i).index  = i;
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

        double getVal(int _index){
            return Entrys.at(_index).Concentration;
        }

        double getLap(int _index){
            return Entrys.at(_index).Lap;
        }

        double getGrad(int _index){
            return Entrys.at(_index).Grad;
        }

        /*************************************************************/

        void addEntry(int num){
            for (int i = 0; i < num; ++i)
            {
                Entrys.push_back(Def_ConEnt);
            }
            Num_Ent = Entrys.size();
            update_index();
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

        void updateEntry(ELEMENT _elemnt,double _con){
            for(auto &ent : Entrys){
                if(ent.Element == _elemnt){
                    ent.Concentration = _con;
                    return;
                }
            }
            throw std::invalid_argument("No such element");
        }

        void updateEntry(int _index,double _con){
            for(auto &ent : Entrys){
                if(ent.index == _index){
                    ent.Concentration = _con;
                    return;
                }
            }
            throw std::invalid_argument("No such element");
        }

        void updateLap(unsigned _index, double _Lap){
            for(auto &ent : Entrys){
                if(ent.index == _index){
                    ent.Lap = _Lap;
                    return;
                }
            }
            throw std::out_of_range("Not in entry list");
        }

        void updateGrad(unsigned _index, double _Grad){
            for(auto &ent : Entrys){
                if(ent.index == _index){
                    ent.Grad = _Grad;
                    return;
                }
            }
            throw std::out_of_range("Not in entry list");
        }
    
        /*************************************************************/

        void update_index(){
            for(unsigned i = 0; i < Num_Ent; i++){
                Entrys.at(i).index = (i);
            }
        }
}Def_ConNode;

#endif