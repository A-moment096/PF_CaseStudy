#pragma once
#ifndef CON_ENTRY
#define CON_ENTRY

#include <iostream>
#include <vector>

enum class ELEMENT{VA,
H,                                                                  He, 
Li, Be,                                         B,  C,  N,  O,  F,  Ne,
Na, Mg,                                         Al, Si, P,  S,  Cl, Ar, 
K,  Ca, Sc, Ti, V,  Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr,
Rb, Sr, Y,  Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I,  Xe,
Cs, Ba, La, Hf, Ta, W,  Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, 
Fr, Ra, Ac, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Nh, Fl, Mc, Lv, Ts, Og,
        
        Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu,
        Th, Pa, U,  Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr,
};

class ConEntry{
    private:
    double Concentration = 0;
    unsigned index;
    public:
    
    ConEntry(){};
    ~ConEntry(){};
    ConEntry& operator= (const ConEntry& newEntry){
            Concentration = newEntry.Concentration;
            return *this;
        }

    double getCon(){
        return Concentration;
    }
    int getindex(){
        return index;
    }
    void setCon(const double newOrderPara ){
        Concentration = newOrderPara;
    }
    void setindex(const int _index){
        index = _index;
    }

}DefaultConEnt;

#endif