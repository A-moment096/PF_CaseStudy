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

class BaseEntry{
    public:
    int Index = 0;
    double Val = 0;
    double Lap = 0;
    double Grad = 0;

    /*************************************************************/

    BaseEntry():Index(0){}

    BaseEntry(double _val):Index(0), Val(_val){}

    BaseEntry(const BaseEntry &_Entry){
        Val = _Entry.Val;
        Index = _Entry.Index;
        Lap = _Entry.Lap;
        Grad = _Entry.Grad;
    }

    BaseEntry &operator= (const BaseEntry &_Entry){
        Val = _Entry.Val;
        Index = _Entry.Index;
        Lap = _Entry.Lap;
        Grad = _Entry.Grad;
        return *this;
    }

    ~BaseEntry(){
    }


    /**/    bool operator== (const BaseEntry &_Entry){
        bool result;
        ((Val == _Entry.Val)&&(Index == _Entry.Index)) ? (result = true) : (result = false);
        return result;
    }

    /*************************************************************/


};

class ConEntry : public BaseEntry{
    public:
    ELEMENT Element;
    using BaseEntry::BaseEntry;

}Def_ConEnt;

class PhaseEntry : public BaseEntry{
    public:
    using BaseEntry::BaseEntry;

}Def_PhsEnt;


#endif