#pragma once
#ifndef CON_NODE
#define CON_NODE

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include "ConEntry.hh"

class ConNode{
    private:
        unsigned Num_Elemt;
        std::vector<ConEntry> Entrys;
    public:
// Construct & Deconstruct Functions
        ConNode(){
            Entrys.push_back(Def_ConEnt);
            Num_Elemt = 1;
        } //Accept Default Parameters
        
        ConNode(std::vector<ConEntry> EntryList){
            Entrys = EntryList;
            Num_Elemt = EntryList.size();
        }

        ConNode(std::vector<ELEMENT> ElementList){
            Num_Elemt = ElementList.size();
            for (unsigned i = 0; i < Num_Elemt; i++){
                Entrys.push_back(Def_ConEnt);
                Entrys.at(i).setElement(ElementList.at(i));
            }
        }

        ConNode(const ConNode& _node){
            Num_Elemt = _node.Num_Elemt;
            Entrys = _node.Entrys;
        }

        ConNode& operator= (const ConNode& _node){
            Num_Elemt = _node.Num_Elemt;
            Entrys = _node.Entrys;
            return *this;
        }

        ~ConNode(){
            Entrys.clear();
        }

        unsigned getNums(){
            return Num_Elemt;
        }

        ELEMENT getElement(unsigned index){
            return Entrys.at(index).getElement();
        }

        std::vector<ELEMENT> getElementList(){
            std::vector<ELEMENT> result;
            for(auto ent : Entrys){
                result.push_back(ent.getElement());
            }
            return result;
        }

        std::vector<double> getCon(){
            std::vector<double> result;
            for(auto ent : Entrys){
                result.push_back(ent.getCon());
            }
            return result;
        }

        double getCon(ELEMENT _elemnt){
            for(auto ent : Entrys){
                if(ent.getElement()==_elemnt){
                    return ent.getCon();
                }
            }
            throw std::invalid_argument("No such element");
            return 0.0;
        }

        ConEntry getEntry(ELEMENT _elemnt){
            for(auto ent : Entrys){
                if(ent.getElement()==_elemnt){
                    return ent;
                }
            }
            throw std::invalid_argument("No such element");
        }

        void addEntry(ConEntry _entry){
            Entrys.push_back(_entry);
        }

        void addEntry(ELEMENT _element){
            ConEntry tempEntry;
            tempEntry.setElement(_element);
            Entrys.push_back(tempEntry);
        }

        void addEntry(ELEMENT _element, double _con){
            ConEntry tempEntry;
            tempEntry.setElement(_element);
            tempEntry.setCon(_con);
            Entrys.push_back(tempEntry);
        }

        void deletEntry(ELEMENT delElemnt){
            if(Num_Elemt <= 1) throw std::out_of_range("Last entry");
            else
                for(auto &ent : Entrys){
                    if(ent.getElement() == delElemnt){
                        std::vector<ConEntry>::iterator delPos = std::find(Entrys.begin(),Entrys.end(),ent);
                        Entrys.erase(delPos);
                        return;
                    }
                }
            throw std::invalid_argument("No such element");
        }

        void updateEntry(ELEMENT _elemnt,double _con){
            for(auto &ent : Entrys){
                if(ent.getElement() == _elemnt){
                    ent.setCon(_con);
                    return;
                }
            }
            throw std::invalid_argument("No such element");
        }

        void updateEntry(double _con){
            if(Num_Elemt == 1){
                Entrys.at(0).setCon(_con);
                return;
            }
            else throw std::invalid_argument("Not only one element");
        }
}Def_ConNode;

#endif