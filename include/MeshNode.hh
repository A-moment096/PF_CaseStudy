#pragma once
#ifndef MESH_NODE_HH
#define MESH_NODE_HH

#include <iomanip>
#include <iostream>
#include <vector>

#include "BaseNode.hh"
/***********************************************
This is MeshNode Class, A node manager.
Every property-node is controlled here.
************************************************/

enum WHICHPARA { CON,
                 PHSFRAC,
                 TEMP,
                 CUSTOM = 99 };
enum WHICHDIR { DirF,
                DirB,
                DirL,
                DirR,
                DirU,
                DirD };

class MeshNode {
private:
    MeshNode *Up = nullptr;
    MeshNode *Forward = nullptr;
    MeshNode *Backward = nullptr;
    MeshNode *Left = nullptr;
    MeshNode *Right = nullptr;
    MeshNode *Down = nullptr;

public:
    PhaseNode Phs_Node;
    ConNode Con_Node;
    CustNode Cust_Node;
    TempNode Temp_Node;

    /*************************************************************/
    // Construct & Deconstruct Functions

    MeshNode() {
        Phs_Node = Def_PhsNode;
        Con_Node = Def_ConNode;
        Cust_Node = Def_CustNode;
        Temp_Node = Def_TempNode;
    };

    MeshNode(PhaseNode _phs_node, ConNode _con_node) {
        Phs_Node = _phs_node;
        Con_Node = _con_node;
        Cust_Node = Def_CustNode;
        Temp_Node = Def_TempNode;
    }

    MeshNode(PhaseNode _phs_node) : MeshNode(_phs_node, Def_ConNode) {}

    MeshNode(ConNode _con_node) : MeshNode(Def_PhsNode, _con_node) {}

    ~MeshNode() {
        Up = nullptr;
        Down = nullptr;
        Forward = nullptr;
        Backward = nullptr;
        Left = nullptr;
        Right = nullptr;
    };

    /*************************************************************/
    // Manipulate Methods

    void setNbhd(WHICHDIR whichdir, MeshNode *_node) {
        switch (whichdir) {
        case DirF:
            Forward = _node;
            break;
        case DirB:
            Backward = _node;
            break;
        case DirL:
            Left = _node;
            break;
        case DirR:
            Right = _node;
            break;
        case DirU:
            Up = _node;
            break;
        case DirD:
            Down = _node;
            break;
        default:
            break;
        }
    }

    MeshNode *getNbhd(WHICHDIR whichdir) {
        switch (whichdir) {
        case DirF:
            return Forward;
            break;
        case DirB:
            return Backward;
            break;
        case DirL:
            return Left;
            break;
        case DirR:
            return Right;
            break;
        case DirU:
            return Up;
            break;
        case DirD:
            return Down;
            break;
        default:
            throw std::invalid_argument("No such Para");
            return nullptr;
            break;
        }
    }

    /*************************************************************/
    void addEnt(WHICHPARA whichpara, int _num) {
        switch (whichpara) {
        case WHICHPARA::CON:
            Con_Node.addEntry(_num);
            return;
            break;
        case WHICHPARA::PHSFRAC:
            Phs_Node.addEntry(_num);
            return;
            break;
        case WHICHPARA::TEMP:
            Temp_Node.addEntry(_num);
            return;
            break;
        case WHICHPARA::CUSTOM:
            Cust_Node.addEntry(_num);
            return;
            break;
        default:
            throw std::invalid_argument("No such Para");
            return;
            break;
        }
    }

    int getNum_Ent(WHICHPARA whichpara) {
        switch (whichpara) {
        case WHICHPARA::CON:
            return Con_Node.Num_Ent;
            break;
        case WHICHPARA::PHSFRAC:
            return Phs_Node.Num_Ent;
            break;
        case WHICHPARA::TEMP:
            return Temp_Node.Num_Ent;
            break;
        case WHICHPARA::CUSTOM:
            return Cust_Node.Num_Ent;
            break;
        default:
            throw std::invalid_argument("No such Para");
            return 1;
            break;
        }
    }

    double getVal(WHICHPARA whichpara, int _Index) {
        switch (whichpara) {
        case WHICHPARA::CON:
            return Con_Node.getVal(_Index);
            break;
        case WHICHPARA::PHSFRAC:
            return Phs_Node.getVal(_Index);
            break;
        case WHICHPARA::TEMP:
            return Temp_Node.getVal(_Index);
            break;
        case WHICHPARA::CUSTOM:
            return Cust_Node.getVal(_Index);
            break;
        default:
            throw std::invalid_argument("No such Para");
            return 1;
            break;
        }
        return 1;
    }

    double getWeight(WHICHPARA whichpara, int _Index) {
        switch (whichpara) {
        case WHICHPARA::CON:
            return Con_Node.getWeight(_Index);
            break;
        case WHICHPARA::PHSFRAC:
            return Phs_Node.getWeight(_Index);
            break;
        case WHICHPARA::TEMP:
            return Temp_Node.getWeight(_Index);
            break;
        case WHICHPARA::CUSTOM:
            return Cust_Node.getWeight(_Index);
            break;
        default:
            throw std::invalid_argument("No such Para");
            return 1;
            break;
        }
        return 1;
    }

    double getLap(int whichpara, int _Index) {
        switch (whichpara) {
        case WHICHPARA::CON:
            return Con_Node.getLap(_Index);
            break;
        case WHICHPARA::PHSFRAC:
            return Phs_Node.getLap(_Index);
            break;
        case WHICHPARA::TEMP:
            return Temp_Node.getLap(_Index);
            break;
        case WHICHPARA::CUSTOM:
            return Cust_Node.getLap(_Index);
            break;
        default:
            throw std::invalid_argument("No such Para or Index");
            return 1;
            break;
        }
        return 1;
    }

    std::vector<double> getGrad(int whichpara, int _Index) {
        switch (whichpara) {
        case WHICHPARA::CON:
            return Con_Node.getGrad(_Index);
            break;
        case WHICHPARA::PHSFRAC:
            return Phs_Node.getGrad(_Index);
            break;
        case WHICHPARA::TEMP:
            return Temp_Node.getGrad(_Index);
            break;
        case WHICHPARA::CUSTOM:
            return Cust_Node.getGrad(_Index);
            break;
        default:
            throw std::invalid_argument("No such Para or Index");
            return {};
            break;
        }
        return {};
    }

    double getGrad(int whichpara, int _Index, DIM whichdim) {
        switch (whichpara) {
        case WHICHPARA::CON:
            return Con_Node.getGrad(_Index, whichdim);
            break;
        case WHICHPARA::PHSFRAC:
            return Phs_Node.getGrad(_Index, whichdim);
            break;
        case WHICHPARA::TEMP:
            return Temp_Node.getGrad(_Index, whichdim);
            break;
        case WHICHPARA::CUSTOM:
            return Cust_Node.getGrad(_Index, whichdim);
            break;
        default:
            throw std::invalid_argument("No such Para or Index");
            return 1;
            break;
        }
        return {};
    }

    std::vector<double> getVelo(int whichpara, int _Index) {
        switch (whichpara) {
        case WHICHPARA::CON:
            return Con_Node.getVelo(_Index);
            break;
        case WHICHPARA::PHSFRAC:
            return Phs_Node.getVelo(_Index);
            break;
        case WHICHPARA::TEMP:
            return Temp_Node.getVelo(_Index);
            break;
        case WHICHPARA::CUSTOM:
            return Cust_Node.getVelo(_Index);
            break;
        default:
            throw std::invalid_argument("No such Para or Index");
            return {};
            break;
        }
        return {};
    }

    double getVelo(int whichpara, int _Index, DIM whichdim) {
        switch (whichpara) {
        case WHICHPARA::CON:
            return Con_Node.getVelo(_Index, whichdim);
            break;
        case WHICHPARA::PHSFRAC:
            return Phs_Node.getVelo(_Index, whichdim);
            break;
        case WHICHPARA::TEMP:
            return Temp_Node.getVelo(_Index, whichdim);
            break;
        case WHICHPARA::CUSTOM:
            return Cust_Node.getVelo(_Index, whichdim);
            break;
        default:
            throw std::invalid_argument("No such Para or Index");
            return 1;
            break;
        }
        return {};
    }

    void iterateVal(WHICHPARA whichpara, double dtime) {
        switch (whichpara) {
        case WHICHPARA::CON:
            Con_Node.iterateVal(dtime);
            return;
            break;
        case WHICHPARA::PHSFRAC:
            Phs_Node.iterateVal(dtime);
            return;
            break;
        case WHICHPARA::TEMP:
            Temp_Node.iterateVal(dtime);
            return;
            break;
        case WHICHPARA::CUSTOM:
            Cust_Node.iterateVal(dtime);
            return;
            break;
        default:
            throw std::invalid_argument("No such Para or Index");
            return;
            break;
        }
        return;
    }
    /*************************************************************/

    double sumPhsFrac() {
        return Phs_Node.sumPhsFrac();
    }

    double sumPhsFrac2() {
        return Phs_Node.sumPhsFrac2();
    }

    double sumPhsFrac3() {
        return Phs_Node.sumPhsFrac3();
    }

    /*************************************************************/

    void showNode();

} Def_Node;

/*************************************************************/
void MeshNode::showNode() {
    std::cout << "Node Information:\n";
    std::cout << "Phase Index:\tPhase Fraction:\t\tElement:\tConcentration:\n";
    for (int i = 0; i < Phs_Node.Num_Ent; i++) {
        for (int j = 0; j < Con_Node.Num_Ent; j++) {
            std::cout << Phs_Node.Entrys.at(i).Index << "\t\t" << std::fixed << std::setprecision(6) << Phs_Node.getVal(i) << "\t\t";
            std::cout << Con_Node.Entrys.at(j).Element << "\t\t" << std::fixed << std::setprecision(6) << Con_Node.getVal(j) << "\n";
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
}

#endif
