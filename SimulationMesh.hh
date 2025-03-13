#pragma once
#ifndef SIMULATION_MESH_HH
#define SIMULATION_MESH_HH

#include <iostream>
#include <vector>
#include <fstream>
#include <filesystem>
#include <string>
#include <cstring>
#include "MeshNode.hh"
#include "PFMTools.hh"

enum BOUNDCOND{ BndPERIODIC, BndCONST, BndADIABATIC };
enum STENCILE{ StencFIVE = 5, StencNINE = 9 };


class SimulationMesh{
    public:
    std::vector<int> Dimension{ 64, 64, 1 };
    int &MeshX = Dimension.at(0);
    int &MeshY = Dimension.at(1);
    int &MeshZ = Dimension.at(2);

    std::vector<double> StepLength{ 1, 1, 1 };
    double &StepX = StepLength.at(0);
    double &StepY = StepLength.at(1);
    double &StepZ = StepLength.at(2);

    double TimeStep ;

    int Num_Nodes;

    BOUNDCOND Bound_Cond;

    std::vector<MeshNode> SimuNodes;


    /*************************************************************/

    SimulationMesh() = delete;

    SimulationMesh(std::vector<int> SizeInfo, std::vector<double> StepInfo, double _TimeStep, MeshNode Node){ // initial with size and SimuNodes.at(i)s
        Dimension = SizeInfo;
        StepLength = StepInfo;
        TimeStep = _TimeStep;
        fillNodes(Node);
        Num_Nodes = SimuNodes.size();
        Bound_Cond = BOUNDCOND::BndPERIODIC;
        bindNbhd();
        bindBoundary(Bound_Cond);
    }

    // SimulationMesh(std::vector<double> StepInfo, MeshNode Node, double _TimeStep):SimulationMesh({ 64, 64, 1 }, StepInfo, _TimeStep , Node){}
    // SimulationMesh(MeshNode Node):SimulationMesh({ 1, 1, 1 },  ,Node){} // initial with SimuNodes.at(i)s and default size

    ~SimulationMesh(){
        Dimension.clear();
        StepLength.clear();
        SimuNodes.clear();
    };

    MeshNode &operator()(const int where){

        if (!(where<Num_Nodes)){
            throw std::invalid_argument("_Index Not in Mesh");
        }
        return SimuNodes.at(where);
    }

    MeshNode &operator()(int X, int Y, int Z){
        if (!(X<Dimension.at(0)&&Y<Dimension.at(1)&&Z<Dimension.at(2) && !(X<0) &&!(Y<0) &&!(Z<0))){
            throw std::invalid_argument("_Index Not in Mesh");
        }
        return SimuNodes.at(X+Y*Dimension.at(0)+Z*Dimension.at(0)*Dimension.at(1));
    }

    /*************************************************************/
    //Initialize tools

    void fillNodes(MeshNode Nodes){
        SimuNodes.reserve(MeshX*MeshY*MeshZ);
        for (int i = 0; i < MeshX*MeshY*MeshZ; i++){
            SimuNodes.push_back(Nodes);
        }
    }

    void bindNbhd(){
        for (int i = 0; i < Num_Nodes; i++){
            SimuNodes.at(i).setNbhd(WHICHDIR::DirF, (i-MeshY < 0                ? &(SimuNodes.at(i)) : &(SimuNodes.at(i-MeshY))));
            SimuNodes.at(i).setNbhd(WHICHDIR::DirB, (i+MeshY >= Num_Nodes       ? &(SimuNodes.at(i)) : &(SimuNodes.at(i+MeshY))));
            SimuNodes.at(i).setNbhd(WHICHDIR::DirL, (i-1 < 0                    ? &(SimuNodes.at(i)) : &(SimuNodes.at(i-1))));
            SimuNodes.at(i).setNbhd(WHICHDIR::DirR, (i+1 >= Num_Nodes           ? &(SimuNodes.at(i)) : &(SimuNodes.at(i+1))));
            SimuNodes.at(i).setNbhd(WHICHDIR::DirU, (i-MeshX*MeshY < 0          ? &(SimuNodes.at(i)) : &(SimuNodes.at(i-MeshX*MeshY))));
            SimuNodes.at(i).setNbhd(WHICHDIR::DirD, (i+MeshX*MeshY >= Num_Nodes ? &(SimuNodes.at(i)) : &(SimuNodes.at(i+MeshX*MeshY))));
        }
    }

    void bindBoundary(BOUNDCOND whichBOUNDCOND){
        for (int i = 0; i < Num_Nodes; i++){
            if (whichBOUNDCOND == BOUNDCOND::BndPERIODIC){
                //forward
                if ((i%(MeshX*MeshY))>=0 && (i%(MeshX*MeshY))<MeshX){
                    SimuNodes.at(i).setNbhd(WHICHDIR::DirF, &(SimuNodes.at(i+MeshX*(MeshY-1))));
                }
                //backward
                if ((i%(MeshX*MeshY))>=MeshX*(MeshY-1) && (i%(MeshX*MeshY))<MeshX*MeshY){
                    SimuNodes.at(i).setNbhd(WHICHDIR::DirB, &(SimuNodes.at(i-MeshX*(MeshY-1))));
                }
                //left
                if (i%MeshX == 0){
                    SimuNodes.at(i).setNbhd(WHICHDIR::DirL, &(SimuNodes.at(i+(MeshX-1))));
                }
                //right
                if (i%MeshX == (MeshX-1)){
                    SimuNodes.at(i).setNbhd(WHICHDIR::DirR, &(SimuNodes.at(i-(MeshX-1))));
                }
                //up
                if (i>=0 && i<MeshX*MeshY){
                    SimuNodes.at(i).setNbhd(WHICHDIR::DirU, &(SimuNodes.at(i+MeshX*MeshY*(MeshZ-1))));
                }
                //down
                if (i>=MeshX*MeshY*(MeshZ-1) && i<MeshX*MeshY*MeshZ){
                    SimuNodes.at(i).setNbhd(WHICHDIR::DirD, &(SimuNodes.at(i-MeshX*MeshY*(MeshZ-1))));
                }
            }

            if (whichBOUNDCOND == BOUNDCOND::BndCONST || whichBOUNDCOND == BOUNDCOND::BndADIABATIC){
                //forward
                if ((i%(MeshX*MeshY))>=0 && (i%(MeshX*MeshY))<MeshX){
                    SimuNodes.at(i).setNbhd(WHICHDIR::DirF, &SimuNodes.at(i));
                }
                //backward
                if ((i%(MeshX*MeshY))>=MeshX*(MeshY-1) && (i%(MeshX*MeshY))<MeshX*MeshY){
                    SimuNodes.at(i).setNbhd(WHICHDIR::DirB, &SimuNodes.at(i));
                }
                //left
                if (i%MeshX == 0){
                    SimuNodes.at(i).setNbhd(WHICHDIR::DirL, &SimuNodes.at(i));
                }
                //right
                if (i%MeshX == (MeshX-1)){
                    SimuNodes.at(i).setNbhd(WHICHDIR::DirR, &SimuNodes.at(i));
                }
                //up
                if (i>=0 && i<MeshX*MeshY){
                    SimuNodes.at(i).setNbhd(WHICHDIR::DirU, &SimuNodes.at(i));
                }
                //down
                if (i>=MeshX*MeshY*(MeshZ-1) && i<MeshX*MeshY*MeshZ){
                    SimuNodes.at(i).setNbhd(WHICHDIR::DirD, &SimuNodes.at(i));
                }
            }
        }
    }

    void ConstBondInit(WHICHPARA whichpara,std::vector<double> _val){
        if(Bound_Cond == BOUNDCOND::BndCONST){
            for(int num = 0; num < getNum_Ent(whichpara); num++){
                for (int i = 0; i < Num_Nodes; i++){
                    //forward
                    if ((i%(MeshX*MeshY))>=0 && (i%(MeshX*MeshY))<MeshX){
                        updateNodeVal(whichpara,i,num,_val.at(0));
                    }
                    //backward
                    if ((i%(MeshX*MeshY))>=MeshX*(MeshY-1) && (i%(MeshX*MeshY))<MeshX*MeshY){
                        updateNodeVal(whichpara,i,num,_val.at(1));
                    }
                    //left
                    if (i%MeshX == 0){
                        updateNodeVal(whichpara,i,num,_val.at(2));
                    }
                    //right
                    if (i%MeshX == (MeshX-1)){
                        updateNodeVal(whichpara,i,num,_val.at(3));
                    }
                    //up
                    if (i>=0 && i<MeshX*MeshY){
                        updateNodeVal(whichpara,i,num,_val.at(4));
                    }
                    //down
                    if (i>=MeshX*MeshY*(MeshZ-1) && i<MeshX*MeshY*MeshZ){
                        updateNodeVal(whichpara,i,num,_val.at(5));
                    }
                }
            }
        }
    }

    /*************************************************************/

    void addEntry(WHICHPARA whichpara, int num){
        switch (whichpara){
        case WHICHPARA::CON:
            for (auto &node : SimuNodes)node.Con_Node.addEntry(num);
            break;
        case WHICHPARA::PHSFRAC:
            for (auto &node : SimuNodes)node.Phs_Node.addEntry(num);
            break;
        case WHICHPARA::TEMP:
            for (auto &node : SimuNodes)node.Temp_Node.addEntry(num);
            break;
        case WHICHPARA::CUSTOM:
            for (auto &node : SimuNodes)node.Cust_Node.addEntry(num);
            break;
        default:
            break;
        }
    }

    void deletEntry(WHICHPARA whichpara, int _Index){
        switch (whichpara)
        {
        case WHICHPARA::CON:
            for(auto &node : SimuNodes)node.Con_Node.deletEntry(_Index);
            break;
        case WHICHPARA::PHSFRAC:
            for(auto &node : SimuNodes)node.Phs_Node.deletEntry(_Index);
            break;
        case WHICHPARA::TEMP:
            for(auto &node : SimuNodes)node.Temp_Node.deletEntry(_Index);
            break;
        case WHICHPARA::CUSTOM:
            for(auto &node : SimuNodes)node.Cust_Node.deletEntry(_Index);
            break;                    
        default:
            break;
        }
    }
    
    void setTimeStep(double _TimeStep){
        TimeStep = _TimeStep;
    }

    /*************************************************************/

    // get value of whichpara at _Index's entry
    std::vector<double> getMeshProp(WHICHPARA whichpara, int _Index){
        std::vector<double> result(Num_Nodes, 0.0);
        if (getNum_Ent(whichpara) != 0 &&_Index<getNum_Ent(whichpara))
        #pragma omp parallel for
            for (int i = 0; i < Num_Nodes;i++){
                result.at(i) = (SimuNodes.at(i).getVal(whichpara,_Index));
            }
        else throw std::out_of_range("_Index out of range");
        return result;
    }

    std::vector<double> getUni_Prop(WHICHPARA whichpara){
        std::vector<double> result;
        for (auto node : SimuNodes){
            double sum = 0;
            for (int i = 0; i < node.getNum_Ent(whichpara); i++){
                sum += node.getVal(whichpara,i)* node.getVal(whichpara,i);
            }
            PFMTools::threshold(sum);
            result.push_back(sum);
        }
        return result;
    }

    /*************************************************************/
    

    std::vector<int> transCoord(int where){
        if (where<Num_Nodes){
            std::vector<int>coord(3);
            coord.at(0) = where%MeshX;
            coord.at(1) = ((where-coord.at(0)) / MeshX)%MeshY;
            coord.at(2) = (where-(coord.at(0)+coord.at(1)*MeshX)/(MeshX*MeshY));
            return coord;
        }
        throw std::out_of_range("Not in mesh");
        return {};
    }

    int transCoord(std::vector<int> where){
        if (where.at(0)<MeshX&&where.at(1)<MeshY&&where.at(2)<MeshZ)
            return where.at(0)+where.at(1)*MeshX+where.at(2)*MeshX*MeshY;
        throw std::out_of_range("Not in mesh");
        return 0;
    }

    int getNum_Ent(WHICHPARA which){
        return SimuNodes.at(0).getNum_Ent(which);
    }

    /*************************************************************/
    //! These methods is 2D only!//

    bool generateDisk(WHICHPARA whichpara,std::vector<int> coord, int index, double Radius){
        bool flag = true;
        const int minX=(((coord.at(0)-Radius)<0)?(flag = false):(coord.at(0)-Radius));
        const int minY=(((coord.at(1)-Radius)<0)?(flag = false):(coord.at(1)-Radius));
        const int maxX=(((coord.at(0)+Radius)>=MeshX)?(flag = false):(coord.at(0)+Radius));
        const int maxY=(((coord.at(1)+Radius)>=MeshY)?(flag = false):(coord.at(1)+Radius));
        if(flag == false){return flag;}
        for(int i = minX; i <maxX; i++){
            for(int j = minY; j<maxY; j++){
                const int a = coord.at(0)+rand()%3-1;
                const int b = coord.at(0)+rand()%3-1;
                const int c = coord.at(1)+rand()%3-1;
                const int d = coord.at(1)+rand()%3-1;
                if((i-a)*(i-b)+(j-c)*(j-d)< Radius*Radius){
                    updateNodeVal(whichpara,{i,j,0},index,1);
                }
            }
        }
        std::cout<<"one disk generated, index: "<<index<<"\n";
        return flag;
    }

    std::vector<std::vector<int>> gnrtDiskSeeds(int _num, double Radius, double tolerance){
        int generatedNum = 0;
        std::vector<std::vector<int>> result(_num,{0,0});
        while(true){
            bool flag = true;
            const int newX = Radius+rand()%MeshX; 
            const int newY = Radius+rand()%MeshY;
            const int maxX=(((newX+Radius)>=MeshX)?(flag = false):(newX+Radius));
            const int maxY=(((newY+Radius)>=MeshY)?(flag = false):(newY+Radius));
            if(flag == false){break;}
            for(int i  = 0; i < generatedNum; i++){
                // if(result.at(i))
            }
        }
    }
    // bool isOverlap(WHICHPARA whichpara,std::vector<int> newCoord, double Radius, double tolerance){
    //     int flag=0;
    //     std::vector<double> testedVal (getUni_Prop(whichpara));
    //     for(int i = (newCoord.at(0)-Radius-1<0?0:newCoord.at(0)-Radius-1); i < (newCoord.at(0)+Radius+1>=MeshX?MeshX:newCoord.at(0)+Radius+1); i++){
    //         for(int j = (newCoord.at(1)-Radius-1<0?0:newCoord.at(1)-Radius-1); j<(newCoord.at(1)+Radius+1>=MeshY?MeshY:newCoord.at(1)+Radius+1); j++){
    //             if((i-newCoord.at(0))*(i-newCoord.at(0))+(j-newCoord.at(1))*(j-newCoord.at(1))<= (Radius-tolerance)*(Radius-tolerance)){
    //                 if(testedVal.at(transCoord({i,j,0}))>0.00001){flag++;}
    //             }
    //         }
    //     }
    //     testedVal.clear();
    //     if(flag != 0){
    //         return true;
    //     }
    //     else{
    //         return false;
    //     }
    // }

    /*************************************************************/

    void updateNodeVal(WHICHPARA whichpara, int where, int _Index, double _val){
        switch (whichpara)
        {
        case WHICHPARA::CON:
            SimuNodes.at(where).Con_Node.updateVal(_Index,_val);
            break;
        case WHICHPARA::PHSFRAC:
            SimuNodes.at(where).Phs_Node.updateVal(_Index,_val);
            break;
        case WHICHPARA::TEMP:
            SimuNodes.at(where).Temp_Node.updateVal(_Index,_val);
            break;
        case WHICHPARA::CUSTOM:
            SimuNodes.at(where).Cust_Node.updateVal(_Index,_val);
            break;
        default:
            throw std::invalid_argument("No such para");
            break;
        }
    }

    void updateNodeVal(WHICHPARA whichpara, std::vector<int> where, int _Index, double _val){
        updateNodeVal(whichpara,transCoord(where), _Index, _val);
    }

    void iterateVal(WHICHPARA whichpara){
        for(auto &node : SimuNodes){
            node.iterateVal(whichpara,TimeStep);
        }
    }
    /*************************************************************/

    /**/
    void Laplacian(STENCILE whichSTNCL, WHICHPARA whichpara){
        double result = 0;
        int Num_Ent = SimuNodes.at(0).getNum_Ent(whichpara);
        if (whichSTNCL == STENCILE::StencFIVE){
        #pragma omp parallel for collapse(2)
        for (auto &node : SimuNodes){
            for (int _Index = 0; _Index < Num_Ent; ++_Index){
                    double c = node.getVal(whichpara,_Index);
                    double f = node.getNbhd(WHICHDIR::DirF)->getVal(whichpara,_Index);
                    double b = node.getNbhd(WHICHDIR::DirB)->getVal(whichpara,_Index);
                    double l = node.getNbhd(WHICHDIR::DirL)->getVal(whichpara,_Index);
                    double r = node.getNbhd(WHICHDIR::DirR)->getVal(whichpara,_Index);
                    // double u = node.getNbhd(WHICHDIR::DirU)->getVal(whichpara,_Index);
                    // double d = node.getNbhd(WHICHDIR::DirD)->getVal(whichpara,_Index);

                    result = ((f+b+l+r-4*c)/(StepX*StepY*StepZ));
                    switch (whichpara)
                    {
                    case WHICHPARA::CON:
                        node.Con_Node.updateLap(_Index,result);
                        break;
                    case WHICHPARA::PHSFRAC:
                        node.Phs_Node.updateLap(_Index,result);
                        break;
                    case WHICHPARA::TEMP:
                        node.Temp_Node.updateLap(_Index,result);
                        break;
                    case WHICHPARA::CUSTOM:
                        node.Cust_Node. updateLap(_Index,result);
                        break;
                    default:
                        break;
                    }
                }
            }
        }
        return;
    }

    void Laplacian(WHICHPARA whichpara){
        Laplacian(STENCILE::StencFIVE, whichpara);
        return;
    }

    void Gradient(WHICHPARA whichpara){
        int Num_Ent = SimuNodes.at(0).getNum_Ent(whichpara);
        #pragma omp parallel for collapse(2)
        for (auto &node : SimuNodes){
            for (int _Index = 0; _Index < Num_Ent; ++_Index){

                    double resultx = (node.getNbhd(WHICHDIR::DirR)->getVal(whichpara,_Index)-node.getNbhd(WHICHDIR::DirL)->getVal(whichpara,_Index))/StepX;
                    double resulty = (node.getNbhd(WHICHDIR::DirB)->getVal(whichpara,_Index)-node.getNbhd(WHICHDIR::DirF)->getVal(whichpara,_Index))/StepY;
                    double resultz = (node.getNbhd(WHICHDIR::DirD)->getVal(whichpara,_Index)-node.getNbhd(WHICHDIR::DirU)->getVal(whichpara,_Index))/StepZ;

                    switch (whichpara)
                    {
                    case WHICHPARA::CON:
                        node.Con_Node.updateGrad(DIM::DimX,_Index,resultx);
                        node.Con_Node.updateGrad(DIM::DimY,_Index,resulty);
                        node.Con_Node.updateGrad(DIM::DimZ,_Index,resultz);
                        break;
                    case WHICHPARA::PHSFRAC:
                        node.Phs_Node.updateGrad(DIM::DimX,_Index,resultx);
                        node.Phs_Node.updateGrad(DIM::DimY,_Index,resulty);
                        node.Phs_Node.updateGrad(DIM::DimZ,_Index,resultz);
                        break;
                    case WHICHPARA::TEMP:
                        node.Temp_Node.updateGrad(DIM::DimX,_Index,resultx);
                        node.Temp_Node.updateGrad(DIM::DimY,_Index,resulty);
                        node.Temp_Node.updateGrad(DIM::DimZ,_Index,resultz);
                        break;
                    case WHICHPARA::CUSTOM:
                        node.Cust_Node.updateGrad(DIM::DimX,_Index,resultx);
                        node.Cust_Node.updateGrad(DIM::DimY,_Index,resulty);
                        node.Cust_Node.updateGrad(DIM::DimZ,_Index,resultz);
                        break;
                    default:
                        break;
                }
            }
        }
    }

    void GradientX(WHICHPARA whichpara, int _Index){
        #pragma omp parallel for
        for(auto &node:SimuNodes){
            double resultx = (node.getNbhd(WHICHDIR::DirR)->getVal(whichpara,_Index)-node.getNbhd(WHICHDIR::DirL)->getVal(whichpara,_Index))/StepX;
            switch (whichpara){
            case WHICHPARA::CON:
                node.Con_Node.updateGrad(DIM::DimX,_Index,resultx);
                break;
            case WHICHPARA::PHSFRAC:
                node.Phs_Node.updateGrad(DIM::DimX,_Index,resultx);
                break;
            case WHICHPARA::TEMP:
                node.Temp_Node.updateGrad(DIM::DimX,_Index,resultx);
                break;
            case WHICHPARA::CUSTOM:
                node.Cust_Node.updateGrad(DIM::DimX,_Index,resultx);
                break;
            default:
                break;
            }
        }
    }

    void GradientY(WHICHPARA whichpara, int _Index){
        #pragma omp parallel for
        for (auto &node : SimuNodes){
            double resulty = (node.getNbhd(WHICHDIR::DirB)->getVal(whichpara,_Index)-node.getNbhd(WHICHDIR::DirF)->getVal(whichpara,_Index))/StepY;
            switch (whichpara){
            case WHICHPARA::CON:
                node.Con_Node.updateGrad(DIM::DimY,_Index,resulty);
                break;
            case WHICHPARA::PHSFRAC:
                node.Phs_Node.updateGrad(DIM::DimY,_Index,resulty);
                break;
            case WHICHPARA::TEMP:
                node.Temp_Node.updateGrad(DIM::DimY,_Index,resulty);
                break;
            case WHICHPARA::CUSTOM:
                node.Cust_Node.updateGrad(DIM::DimY,_Index,resulty);
                break;
            default:
                break;
            }
        }
    }
    void GradientZ(WHICHPARA whichpara, int _Index){
        #pragma omp parallel for
        for (auto &node : SimuNodes){
            double resultz = (node.getNbhd(WHICHDIR::DirD)->getVal(whichpara,_Index)-node.getNbhd(WHICHDIR::DirU)->getVal(whichpara,_Index))/StepZ;
            switch (whichpara){
            case WHICHPARA::CON:
                node.Con_Node.updateGrad(DIM::DimZ,_Index,resultz);
                break;
            case WHICHPARA::PHSFRAC:
                node.Phs_Node.updateGrad(DIM::DimZ,_Index,resultz);
                break;
            case WHICHPARA::TEMP:
                node.Temp_Node.updateGrad(DIM::DimZ,_Index,resultz);
                break;
            case WHICHPARA::CUSTOM:
                node.Cust_Node.updateGrad(DIM::DimZ,_Index,resultz);
                break;
            default:
                break;
            }
        }
    }    

    void GradientX(WHICHPARA whichpara){
        int Num_Ent = SimuNodes.at(0).getNum_Ent(whichpara);
        for (int _Index = 0; _Index <Num_Ent; ++_Index){
                GradientX(whichpara,_Index);
        }
    }
    void GradientY(WHICHPARA whichpara){
        int Num_Ent = SimuNodes.at(0).getNum_Ent(whichpara);
            for (int _Index = 0; _Index < Num_Ent; ++_Index){
                GradientY(whichpara,_Index);
            }
    }
    void GradientZ(WHICHPARA whichpara){
        int Num_Ent = SimuNodes.at(0).getNum_Ent(whichpara);
            for (int _Index = 0; _Index < Num_Ent; ++_Index){
                GradientZ(whichpara,_Index);
            }
    }

    double CalDensity(){
        double sum;
        for(auto node : SimuNodes)
            for(int i = 0; i < node.getNum_Ent(WHICHPARA::PHSFRAC); i++)
                sum += node.Phs_Node.getVal(i);
        return sum/Num_Nodes;
    }
    /*************************************************************/
    void showGlobalInfo(); 
    void showNodesProp(WHICHPARA which, int _Index); 

    void outVTKFilehead(std::string _dirname, int istep);
    void outVTKAve(std::string _dirname, WHICHPARA whichpara, int istep);
    void outVTKAll(std::string _dirname, WHICHPARA whichpara, int istep);

    void outVTK(std::string _dirname, int istep){
        outVTKFilehead(_dirname, istep);
        outVTKAve(_dirname, WHICHPARA::CON, istep);
        outVTKAve(_dirname, WHICHPARA::PHSFRAC, istep);
    }

    void outCSV(std::string _dirname, std::string _filename,int istep,double _val){

        std::string filename(_dirname+"/"+_filename+".csv");
        
        std::ofstream outfile;
        outfile.open(filename, std::ios::app);

        outfile<<istep<<","<<_val<<std::endl;

        outfile.close();
    }

    template <typename T>
    void outCSV(std::string _dirname, std::string _filename,std::vector<T> _FirstVals,std::vector<T> _SecondVals){

        std::string filename(_dirname+"/"+_filename+".csv");
        std::ofstream outfile;
        outfile.open(filename);
        for(int i = 0; i < (int) _FirstVals.size(); i++){
            outfile<<_FirstVals.at(i)<<","<<_SecondVals.at(i)<<"\n";
        }
        outfile.flush();
        outfile.close();
    }

};


/**************************************************************************************************************************/
/**************************************************************************************************************************/
/**************************************************************************************************************************/


// show the basic information of the mesh
void SimulationMesh::showGlobalInfo(){
    std::cout<<"SimulationMesh Properties:\n";
    std::cout<<"Mesh Size:\t\t"<<MeshX<<"\u0078"<<MeshY<<"\u0078"<<MeshZ<<"\n";
    std::cout<<"Elements:\t\t";
    for (auto ent : SimuNodes.at(0).Con_Node.Entrys)
        std::cout<<ent.Element<<" ";
    std::cout<<"\nNumber of Phase:\t"<<SimuNodes.at(0).getNum_Ent(WHICHPARA::PHSFRAC);
    std::cout<<"\nNumber of Element:\t"<<SimuNodes.at(0).getNum_Ent(WHICHPARA::CON);
    std::cout<<"\nNumber of CustomEntry:\t"<<SimuNodes.at(0).getNum_Ent(WHICHPARA::CUSTOM);

    std::cout<<"\n-----------------------------------------------------------------\n"<<std::endl;
}

// show one of the properties of the SimuNodes.at(i)s inside the mesh
void SimulationMesh::showNodesProp(WHICHPARA which, int _Index){ //which para, _Index of para
    if (which == WHICHPARA::CON)std::cout<<"Concentration of "<<SimuNodes.at(0).Con_Node.Entrys.at(_Index).Element<<"\n";
    if (which == WHICHPARA::PHSFRAC)std::cout<<"Order Parameter of "<<_Index<<" grain\n";
    if (which == WHICHPARA::CUSTOM)std::cout<<"Custom value \n";
    for (long i = 0; i < Num_Nodes; i++){
        std::cout<<std::fixed<<std::setprecision(10)<<
            SimuNodes.at(i).getVal(which,_Index)
            // SimuNodes.at(i).getLap(WHICHPARA::CON, _Index)
            <<" ";
        if (i%MeshX==MeshX-1)std::cout<<"\n";
        if (i%(MeshX*MeshY) == MeshX*MeshY-1)std::cout<<"\n";
    }
    std::cout<<"-----------------------------------------------------------------\n"<<std::endl;
}

/*
 ! under present folder create a folder 'Result/_dirname' to save vtk files
 ! if "../" is used, the 'Result' will be  omitted
*/
std::string toVTK_Path(std::string _dirname){
    // legal name check
    if (_dirname.empty())throw std::invalid_argument("_dirname can't be empty");
    for (auto _ch : _dirname){
        if (std::isspace(_ch))throw std::invalid_argument("_dirname can't contain space, maybe use underscore");
    }

    std::string result("./Result/");
    result = result + _dirname;
    std::filesystem::create_directories(result);
    return result;
}

inline void SimulationMesh::outVTKFilehead(std::string _dirname, int istep){

    std::ostringstream filename_ss;
    filename_ss<<_dirname<<"/time_"<<istep<<".vtk";
    std::string filename {filename_ss.str()};

    /*************************************************************/
    std::ofstream outfile;
    outfile.open(filename);

    // head
    outfile<<"# vtk DataFile Version 2.0\n";
    outfile<<"time_10.vtk\n";
    outfile<<"ASCII\n";
    outfile<<"DATASET STRUCTURED_GRID\n";

    // size info
    outfile<<"DIMENSIONS "<<MeshX<<"  "<<MeshY<<"  "<<MeshZ<<"\n";
    outfile<<"POINTS "<<Num_Nodes<<"   float\n";
    for (int i = 0; i < MeshX;i++)
        for (int j = 0; j < MeshY;j++)
            for (int k = 0; k < MeshZ; k++){
                outfile<<i*StepX<<"   "<<j*StepY<<"   "<<k*StepZ<<"\n";
            }
    outfile<<"POINT_DATA "<<Num_Nodes<<"\n";

    outfile.flush();
    outfile.close();
}

inline void SimulationMesh::outVTKAve(std::string _dirname, WHICHPARA whichpara, int istep){

    std::ostringstream filename_ss;
    filename_ss<<_dirname<<"/time_"<<istep<<".vtk";
    std::string filename {filename_ss.str()};

    std::ofstream outfile;
    outfile.open(filename, std::ios::app);

    char varname[64];
    switch (whichpara){
    case WHICHPARA::PHSFRAC:
        std::sprintf(varname, "PHSFRAC_AVE");
        break;
    case WHICHPARA::CON:
        std::sprintf(varname, "CON_AVE");
        break;
    case WHICHPARA::CUSTOM:
        std::sprintf(varname, "CUST_AVE");
        break;
    case WHICHPARA::TEMP:
        std::printf(varname, "TEMP_AVE");
        break;
    default:
        break;
    }
    std::vector<double> normed_phs(getUni_Prop(whichpara));
    outfile<<"SCALARS "<<varname<<"  float  1\n";
    outfile<<"LOOKUP_TABLE default\n";
    for (int i = 0;i<Num_Nodes;i++){
        if(!std::isnan(normed_phs.at(i))){
            outfile<<normed_phs.at(i)<<"\n";
        }
        else{
            throw std::invalid_argument("output is nan");
        }
    }
    outfile.flush();

    outfile.close();

}

inline void SimulationMesh::outVTKAll(std::string _dirname, WHICHPARA whichpara, int istep){

    std::ostringstream filename_ss;
    filename_ss<<_dirname<<"/time_"<<istep<<".vtk";
    std::string filename {filename_ss.str()};

    std::ofstream outfile;
    outfile.open(filename, std::ios::app);

    char varname[64];

    for (int num = 0; num < getNum_Ent(whichpara); num++){
        switch (whichpara){
        case WHICHPARA::PHSFRAC:
            std::sprintf(varname, "PHSFRAC_%01d", num);
            break;
        case WHICHPARA::CON:
            std::sprintf(varname, "CON_%01d", num);
            break;
        case WHICHPARA::CUSTOM:
            std::sprintf(varname, "CUST_%01d", num);
            break;
        case WHICHPARA::TEMP:
            std::sprintf(varname, "TEMP_%01d", num);
            break;
        default:
            break;
        }
        outfile<<"SCALARS "<<varname<<"  float  1\n";
        outfile<<"LOOKUP_TABLE default\n";
        for (int i = 0;i<Num_Nodes;i++){
            if(!std::isnan(SimuNodes.at(i).getVal(whichpara,num))){
                outfile<<SimuNodes.at(i).getVal(whichpara,num)<<"\n";
            }
        else{
            throw std::invalid_argument("output is nan");
        }
    }
    }
    outfile.flush();
    outfile.close();
}


#endif