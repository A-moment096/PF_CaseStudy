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
    int Num_Ent;
    std::vector<EntryT> Entrys;

    /*************************************************************/
    
    BaseNode(){
        Entrys.push_back(Def_Ent);
        Num_Ent = 1;
        Entrys.at(0).Index = 0;
    }

    BaseNode(std::vector<EntryT> EntryList){
        Entrys = EntryList;
        Num_Ent = EntryList.size();
        updateIndex(0);
    }

    BaseNode(const BaseNode &_node){
        Num_Ent = _node.Num_Ent;
        Entrys = _node.Entrys;
    }

    BaseNode &operator= (const BaseNode &_node){
        Num_Ent = _node.Num_Ent;
        Entrys = _node.Entrys;
        return *this;
    }

    EntryT &operator() (const int _Index){
        if (_Index<Num_Ent)
            return Entrys.at(_Index);
        else throw std::out_of_range("No such entry");
    }

    ~BaseNode(){
        Entrys.clear();
    }

    /*************************************************************/

    void addEntry(int num){
        for (int i = 0; i < num; ++i){
            Entrys.push_back(Def_Ent);
        }
        Num_Ent = Entrys.size();
        updateIndex(Num_Ent-num);
    }

    void deletEntry(int _Index){
        if (Num_Ent <= 1) throw std::out_of_range("Last entry");
        else if(_Index<Num_Ent && _Index >=0){
            Entrys.erase(Entrys.begin()+_Index);
            updateIndex(_Index);
            return;
        }
        throw std::out_of_range("Not in entry list");
    }

    /*************************************************************/

    std::vector<double> getVal(){
        std::vector<double> result(Num_Ent);
        for (auto ent : Entrys){
            result.push_back(ent.Val);
        }
        return result;
    }

    std::vector<double> getLap(){
        std::vector<double> result(Num_Ent);
        for (auto ent : Entrys){
            result.push_back(ent.Lap);
        }
        return result;
    }

    std::vector<double> getGrad(){
        std::vector<double> result(Num_Ent);
        for (auto ent : Entrys){
            result.push_back(ent.Grad);
        }
        return result;
    }

    double getVal(int _Index){
        if (_Index < Num_Ent)
            return Entrys.at(_Index).Val;
        else throw std::out_of_range("No Such Index");
    }

    double getLap(int _Index){
        if (_Index < Num_Ent)
            return Entrys.at(_Index).Lap;
        else throw std::out_of_range("No Such Index");
    }

    double getGrad(int _Index){
        if (_Index < Num_Ent)
            return Entrys.at(_Index).Grad;
        else throw std::out_of_range("No Such Index");
    }

    /*************************************************************/

    void updateVal(int _Index, double _con){
        if(_Index < Num_Ent && _Index >= 0)
        Entrys.at(_Index).Val = _con;
        else throw std::out_of_range("Not in entry list");
    }

    void updateLap(int _Index, double _Lap){
        if(_Index < Num_Ent && _Index >= 0)
        Entrys.at(_Index).Lap = _Lap;
        else throw std::out_of_range("Not in entry list");
    }

    void updateGrad(int _Index, double _Grad){
        if(_Index < Num_Ent && _Index >= 0)
        Entrys.at(_Index).Grad = _Grad;
        else throw std::out_of_range("Not in entry list");
    }

    /*************************************************************/

    private:
    void updateIndex(int num){
        Num_Ent = Entrys.size();
        for (int i = num; i < Num_Ent; i++){
            Entrys.at(i).Index = i;
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
        for (int i = 0; i < Num_Ent; i++){
            Entrys.at(i).Index = (i);
        }
    }

    double sumPhsFrac(){
        double result = 0;
        for (auto ent : Entrys)result += ent.Val;
        return result;
    }

    double sumPhsFrac2(){
        double result = 0;
        for (auto ent : Entrys)result += ent.Val*ent.Val;
        return result;
    }

    double sumPhsFrac3(){
        double result = 0;
        for (auto ent : Entrys)result += ent.Val*ent.Val*ent.Val;
        return result;
    }
}Def_PhsNode;

class CustNode : public BaseNode<CustEntry>{
    public:
    using BaseNode<CustEntry>::BaseNode;
}Def_CustNode;

class TempNode : public BaseNode<TempEntry>{
    public:
    using BaseNode<TempEntry>::BaseNode;
}Def_TempNode;

#endif
