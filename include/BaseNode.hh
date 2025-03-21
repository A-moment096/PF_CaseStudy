#pragma once
#ifndef NODE_TEMP_HH
#define NODE_TEMP_HH

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

#include "BaseEntry.hh"
#include "PFMTools.hh"

enum DIM { DimX,
           DimY,
           DimZ };

template <class EntryT>
class BaseNode {
private:
    EntryT Def_Ent;

public:
    int Num_Ent;
    std::vector<EntryT> Entrys;

    /*************************************************************/

    BaseNode() {
        Entrys.push_back(Def_Ent);
        Num_Ent = 1;
        Entrys.at(0).Index = 0;
    }

    BaseNode(std::vector<EntryT> EntryList) {
        Entrys = EntryList;
        Num_Ent = EntryList.size();
        updateIndex(0);
    }

    BaseNode(const BaseNode &_node) {
        Num_Ent = _node.Num_Ent;
        Entrys = _node.Entrys;
    }

    BaseNode &operator=(const BaseNode &_node) {
        Num_Ent = _node.Num_Ent;
        Entrys = _node.Entrys;
        return *this;
    }

    EntryT &operator()(const int _Index) {
        if (_Index < Num_Ent)
            return Entrys.at(_Index);
        else
            throw std::out_of_range("No such entry");
    }

    ~BaseNode() {
        Entrys.clear();
    }

    /*************************************************************/

    void addEntry(int num) {
        for (int i = 0; i < num; ++i) {
            Entrys.push_back(Def_Ent);
        }
        Num_Ent = Entrys.size();
        updateIndex(Num_Ent - num);
    }

    void deletEntry(int _Index) {
        if (Num_Ent <= 1)
            throw std::out_of_range("Last entry");
        else if (_Index < Num_Ent && _Index >= 0) {
            Entrys.erase(Entrys.begin() + _Index);
            updateIndex(_Index);
            return;
        }
        throw std::out_of_range("Not in entry list");
    }

    /*************************************************************/

    double getVal(int _Index) {
        if (_Index < Num_Ent) {
            return Entrys.at(_Index).Val;
        } else {
            throw std::out_of_range("No Such Index");
        }
    }

    double getWeight(int _Index) {
        if (_Index < Num_Ent) {
            return Entrys.at(_Index).Weight;
        } else {
            throw std::out_of_range("No Such Index");
        }
    }

    double getLap(int _Index) {
        if (_Index < Num_Ent) {
            return Entrys.at(_Index).Lap;
        } else {
            throw std::out_of_range("No Such Index");
        }
    }

    double getGrad(int _Index, DIM whichdim) {
        if (_Index < Num_Ent) {
            switch (whichdim) {
            case DIM::DimX:
                return Entrys.at(_Index).GradX;
                break;
            case DIM::DimY:
                return Entrys.at(_Index).GradY;
                break;
            case DIM::DimZ:
                return Entrys.at(_Index).GradZ;
                break;
            default:
                break;
            }
        } else {
            throw std::out_of_range("No Such Index");
        }
        return 1;
    }

    std::vector<double> getGrad(int _Index) {
        if (_Index < Num_Ent) {
            return Entrys.at(_Index).Grad;
        } else {
            throw std::out_of_range("No Such Index");
        }
    }

    double getVelo(int _Index, DIM whichdim) {
        if (_Index < Num_Ent) {
            switch (whichdim) {
            case DIM::DimX:
                return Entrys.at(_Index).VeloX;
                break;
            case DIM::DimY:
                return Entrys.at(_Index).VeloY;
                break;
            case DIM::DimZ:
                return Entrys.at(_Index).VeloZ;
                break;
            default:
                break;
            }
        } else {
            throw std::out_of_range("No Such Index");
        }
        return 1;
    }

    std::vector<double> getVelo(int _Index) {
        if (_Index < Num_Ent) {
            return Entrys.at(_Index).Velo;
        } else {
            throw std::out_of_range("No Such Index");
        }
    }

    double getDVal(int _Index) {
        if (_Index < Num_Ent) {
            return Entrys.at(_Index).DVal;
        } else {
            throw std::out_of_range("No Such Index");
        }
    }
    /*************************************************************/

    void updateVal(int _Index, double _val) {
        if (_Index < Num_Ent && _Index >= 0) {
            Entrys.at(_Index).Val = _val;
        } else {
            throw std::out_of_range("Not in entry list");
        }
    }

    void updateWeight(int _Index, double _weight) {
        if (_Index < Num_Ent && _Index >= 0) {
            Entrys.at(_Index).Weight = _weight;
        } else {
            throw std::out_of_range("Not in entry list");
        }
    }

    void updateLap(int _Index, double _Lap) {
        if (_Index < Num_Ent && _Index >= 0) {
            Entrys.at(_Index).Lap = _Lap;
        } else {
            throw std::out_of_range("Not in entry list");
        }
    }

    void updateGrad(DIM whichdim, int _Index, double _Grad) {
        if (_Index < Num_Ent && _Index >= 0) {
            switch (whichdim) {
            case DIM::DimX:
                Entrys.at(_Index).GradX = _Grad;
                break;
            case DIM::DimY:
                Entrys.at(_Index).GradY = _Grad;
                break;
            case DIM::DimZ:
                Entrys.at(_Index).GradZ = _Grad;
                break;
            default:
                break;
            }
        } else {
            throw std::out_of_range("Not in entry list");
        }
    }

    void updateGrad(int _Index, std::vector<double> _Grad) {
        if (_Index < Num_Ent && _Index >= 0) {
            Entrys.at(_Index).Grad = _Grad;
        } else {
            throw std::out_of_range("Not in entry list");
        }
    }

    void updateVelo(DIM whichdim, int _Index, double _Velo) {
        if (_Index < Num_Ent && _Index >= 0) {
            switch (whichdim) {
            case DIM::DimX:
                Entrys.at(_Index).VeloX = _Velo;
                break;
            case DIM::DimY:
                Entrys.at(_Index).VeloY = _Velo;
                break;
            case DIM::DimZ:
                Entrys.at(_Index).VeloZ = _Velo;
                break;
            default:
                break;
            }
        } else {
            throw std::out_of_range("Not in entry list");
        }
    }

    void updateVelo(int _Index, std::vector<double> _Velo) {
        if (_Index < Num_Ent && _Index >= 0) {
            Entrys.at(_Index).Velo = _Velo;
        } else {
            throw std::out_of_range("Not in entry list");
        }
    }

    void updateDVal(int _Index, double _DVal) {
        if (_Index < Num_Ent && _Index >= 0) {
            Entrys.at(_Index).DVal = _DVal;
        } else {
            throw std::out_of_range("Not in entry list");
        }
    }

    void iterateVal(double dtime) noexcept(true) {
        for (auto &ent : Entrys) {
            ent.Val = PFMTools::threshold(ent.Val + dtime * ent.DVal);
        }
    }

    /*************************************************************/

private:
    void updateIndex(int num) {
        Num_Ent = Entrys.size();
        for (int i = num; i < Num_Ent; i++) {
            Entrys.at(i).Index = i;
        }
    }
};

class ConNode : public BaseNode<ConEntry> {
public:
    using BaseNode<ConEntry>::BaseNode;
} Def_ConNode;

class PhaseNode : public BaseNode<PhaseEntry> {
public:
    using BaseNode<PhaseEntry>::BaseNode;
    void updateIndex() {
        for (int i = 0; i < Num_Ent; i++) {
            Entrys.at(i).Index = (i);
        }
    }

    double sumPhsFrac() {
        double result = 0;
        for (auto ent : Entrys)
            result += ent.Val;
        return result;
    }

    double sumPhsFrac2() {
        double result = 0;
        for (auto ent : Entrys)
            result += ent.Val * ent.Val;
        return result;
    }

    double sumPhsFrac3() {
        double result = 0;
        for (auto ent : Entrys)
            result += ent.Val * ent.Val * ent.Val;
        return result;
    }
} Def_PhsNode;

class CustNode : public BaseNode<CustEntry> {
public:
    using BaseNode<CustEntry>::BaseNode;
} Def_CustNode;

class TempNode : public BaseNode<TempEntry> {
public:
    using BaseNode<TempEntry>::BaseNode;
} Def_TempNode;

#endif
