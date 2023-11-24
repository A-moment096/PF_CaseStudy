#pragma once
#ifndef NODE_TEMP_HH
#define NODE_TEMP_HH

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include "BaseEntry.hh"

template <class EntryT>
class BaseNode{
    private:
        EntryT Def_Ent;
    public:
        unsigned Num_Ent;
        std::vector<EntryT> Entrys;

/*************************************************************/

// Construct & Deconstruct Functions
        BaseNode(){
            Entrys.push_back(Def_Ent);
            Num_Ent = 1;
            (*this)(0).Index  = 0;
        } //Accept Default Parameters
        
        BaseNode(std::vector<EntryT> EntryList){
            Entrys = EntryList;
            Num_Ent = EntryList.size();
            updateIndex();
        }

        BaseNode(const BaseNode& _node){
            Num_Ent = _node.Num_Ent;
            Entrys = _node.Entrys;
        }

        BaseNode& operator= (const BaseNode& _node){
            Num_Ent = _node.Num_Ent;
            Entrys = _node.Entrys;
            return *this;
        }
                
        EntryT& operator() (const unsigned _Index){
            if(_Index<Num_Ent)
            return Entrys.at(_Index);
            else throw std::out_of_range("No such entry");
        }

        ~BaseNode(){
            Entrys.clear();
        }

        /*************************************************************/

        std::vector<double> getVal(){
            std::vector<double> result;
            for(auto ent : Entrys){
                result.push_back(ent.Val);
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
            return (*this)(_Index).Val;
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
                Entrys.push_back(Def_Ent);
            }
            Num_Ent = Entrys.size();
            updateIndex();
        }

        void deletEntry(unsigned _Index){
            if(Num_Ent <= 1) throw std::out_of_range("Last entry");
            else
                for(auto &ent : Entrys){
                    if(ent.Index == _Index){
                        std::iterator delPos = std::find(Entrys.begin(),Entrys.end(),ent);
                        Entrys.erase(delPos);
                        return;
                    }
                }
            throw std::out_of_range("Not in entry list");
        }

        void updateVal(unsigned _Index,double _con){
            for(auto &ent : Entrys){
                if(ent.Index == _Index){
                    ent.Val = _con;
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
    
};

class ConNode : public BaseNode<ConEntry>{
    public:
        using BaseNode<ConEntry>::BaseNode;
}Def_ConNode;

class PhaseNode : public BaseNode<PhaseEntry>{
    public:
        using BaseNode<PhaseEntry>::BaseNode;
        void updateIndex(){
            for(unsigned i = 0; i < Num_Ent; i++){
                (*this)(i).Index = (i);
            }
        }

        double sumPhsFrac(){
            double result = 0;
            for(auto ent : Entrys)result += ent.Val;
            return result;  
        }

        double sumPhsFrac2(){
            double result = 0;
            for(auto ent : Entrys)result += ent.Val*ent.Val;
            return result;
        }

        double sumPhsFrac3(){
            double result = 0;
            for(auto ent : Entrys)result += ent.Val*ent.Val*ent.Val;
            return result;
        }
}Def_PhsNode;

#endif
