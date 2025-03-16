#pragma once
#ifndef BASE_ENTRY
#define BASE_ENTRY

#include <iostream>
#include <vector>

enum class ELEMENT:int{
VA=0,
H,                                                                  He, 
Li, Be,                                         B,  C,  N,  O,  F,  Ne,
Na, Mg,                                         Al, Si, P,  S,  Cl, Ar, 
K,  Ca, Sc, Ti, V,  Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr,
Rb, Sr, Y,  Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I,  Xe,
Cs, Ba, La, Hf, Ta, W,  Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, 
Fr, Ra, Ac, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Nh, Fl, Mc, Lv, Ts, Og,
        
        Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu,
        Th, Pa, U,  Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr,
ELEMENTCount
};

std::ostream& operator<<(std::ostream& out, const ELEMENT _element){
const std::vector<std::string> element_list = {
    "VA",
    "H" ,                                                                                                 "He", 
    "Li", "Be",                                                             "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg",                                                             "Al", "Si", "P",  "S",  "Cl", "Ar", 
    "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
    "Cs", "Ba", "La", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", 
    "Fr", "Ra", "Ac", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",

                "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
                "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",

    };
    if(_element >=ELEMENT::ELEMENTCount || _element < ELEMENT::VA)return out<<"???";
    return out<<element_list[(int)_element];
}

class BaseEntry {
public:
    int Index = 0;
    double Val = 0;
    double Weight = 1;

    double Lap = 0;
    std::vector<double> Grad{0.0, 0.0, 0.0};
    double &GradX = Grad.at(0);
    double &GradY = Grad.at(1);
    double &GradZ = Grad.at(2);
    std::vector<double> Velo{0.0, 0.0, 0.0};
    double &VeloX = Velo.at(0);
    double &VeloY = Velo.at(1);
    double &VeloZ = Velo.at(2);

    double DVal = 0;

    /*************************************************************/

    BaseEntry() : Index(0) {}

    BaseEntry(double _val) : Index(0), Val(_val) {}

    BaseEntry(const BaseEntry &_Entry) {
        Index = _Entry.Index;
        Weight = _Entry.Weight;
        Val = _Entry.Val;

        Lap = _Entry.Lap;
        GradX = _Entry.GradX;
        GradY = _Entry.GradY;
        GradZ = _Entry.GradZ;

        VeloX = _Entry.VeloX;
        VeloY = _Entry.VeloY;
        VeloZ = _Entry.VeloZ;

        DVal = _Entry.DVal;
    }

    BaseEntry &operator=(const BaseEntry &_Entry) {
        Index = _Entry.Index;
        Weight = _Entry.Weight;
        Val = _Entry.Val;

        Lap = _Entry.Lap;
        GradX = _Entry.GradX;
        GradY = _Entry.GradY;
        GradZ = _Entry.GradZ;

        VeloX = _Entry.VeloX;
        VeloY = _Entry.VeloY;
        VeloZ = _Entry.VeloZ;

        DVal = _Entry.DVal;
        return *this;
    }

    ~BaseEntry() {
        Grad.clear();
    }

    /*************************************************************/
};

class ConEntry : public BaseEntry {
public:
    ELEMENT Element;
    using BaseEntry::BaseEntry;

} Def_ConEnt;

class PhaseEntry : public BaseEntry {
public:
    using BaseEntry::BaseEntry;

} Def_PhsEnt;

class CustEntry : public BaseEntry {
public:
    using BaseEntry::BaseEntry;
};

class TempEntry : public BaseEntry {
public:
    using BaseEntry::BaseEntry;
} Def_TempEnt;

#endif