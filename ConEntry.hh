#pragma once
#ifndef CON_ENTRY
#define CON_ENTRY

#include <iostream>
#include <string>
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

class ConEntry{
    private:
        double Concentration = 0;
        ELEMENT Element;
        double Lap = 0;
        double Grad = 0;
    public:
    
        ConEntry(){
            Element = ELEMENT::VA;
        }
        ConEntry(ELEMENT _element){
            Element = _element;
        }
        ConEntry(ELEMENT _element, double _con):ConEntry(_element){
            Concentration = _con;
        }

        ConEntry(const ConEntry& _Entry){
            Concentration = _Entry.Concentration;
            Element = _Entry.Element;
        }

        ConEntry& operator= (const ConEntry& _Entry){
            Concentration = _Entry.Concentration;
            Element = _Entry.Element;
            return *this;
        }

        ~ConEntry(){
        }
        
/*************************************************************/

        bool operator== (const ConEntry& _Entry){
            bool result;
            ((Concentration == _Entry.Concentration)&&(Element == _Entry.Element)) ? (result = true) : (result = false);
            return result;
        }

        double getCon(){
            return Concentration;
        }
        ELEMENT getElement(){
            return Element;
        }
        double getLap(){
            return Lap;
        }
        double getGrad(){
            return Grad;
        }


        void setCon(const double _Con ){
            Concentration = _Con;
        }
        void setElement(const ELEMENT _element){
            Element = _element;
        }
        void setLap(const double _Lap){
            Lap = _Lap;
        }
        void setGrad(const double _Grad){
            Grad = _Grad;
        }

}Def_ConEnt;

#endif