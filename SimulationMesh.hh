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

enum DIM{ DimX, DimY, DimZ };
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

    int Num_Nodes;

    std::vector<MeshNode> SimuNodes;


    /*************************************************************/

    SimulationMesh(){} // initial with default size, without SimuNodes.at(i)s

    SimulationMesh(std::vector<int> SizeInfo, std::vector<double> StepInfo, MeshNode Node){ // initial with size and SimuNodes.at(i)s
        Dimension = SizeInfo;
        StepLength = StepInfo;
        fillNodes(Node);
        Num_Nodes = SimuNodes.size();
        bindBoundary(BOUNDCOND::BndPERIODIC);
    }

    SimulationMesh(std::vector<double> StepInfo, MeshNode Node):SimulationMesh({ 64, 64, 1 }, StepInfo, Node){}
    SimulationMesh(MeshNode Node):SimulationMesh({ 1, 1, 1 }, Node){} // initial with SimuNodes.at(i)s and default size

    ~SimulationMesh(){
        Dimension.clear();
        StepLength.clear();
        SimuNodes.clear();
    };

    MeshNode &operator()(const int where){
        if (where<Num_Nodes){
            return SimuNodes.at(where);
        }
        else throw std::invalid_argument("Index Not in Mesh");
    }

    MeshNode &operator()(int X, int Y, int Z){
        if (X<Dimension.at(0)&&Y<Dimension.at(1)&&Z<Dimension.at(2) && !(X<0) &&!(Y<0) &&!(Z<0)){
            return SimuNodes.at(X+Y*Dimension.at(0)+Z*Dimension.at(0)*Dimension.at(1));
        }
        else throw std::invalid_argument("Index Not in Mesh");
    }

    /*************************************************************/

    void fillNodes(MeshNode Nodes){ // fill mesh with SimuNodes.at(i)s 
        SimuNodes.reserve(MeshX*MeshY*MeshZ);
        for (int i = 0; i < MeshX*MeshY*MeshZ; i++){
            SimuNodes.push_back(Nodes);
        }
    }

    void bindBoundary(BOUNDCOND whichBOUNDCOND){
        for (int i = 0; i < Num_Nodes; i++){
            SimuNodes.at(i).setNbhd(WHICHDIR::DirF, (i-MeshY < 0                ? &(SimuNodes.at(i)) : &(SimuNodes.at(i-MeshY))));
            SimuNodes.at(i).setNbhd(WHICHDIR::DirB, (i+MeshY >= Num_Nodes       ? &(SimuNodes.at(i)) : &(SimuNodes.at(i+MeshY))));
            SimuNodes.at(i).setNbhd(WHICHDIR::DirL, (i-1 < 0                    ? &(SimuNodes.at(i)) : &(SimuNodes.at(i-1))));
            SimuNodes.at(i).setNbhd(WHICHDIR::DirR, (i+1 >= Num_Nodes           ? &(SimuNodes.at(i)) : &(SimuNodes.at(i+1))));
            SimuNodes.at(i).setNbhd(WHICHDIR::DirU, (i-MeshX*MeshY < 0          ? &(SimuNodes.at(i)) : &(SimuNodes.at(i-MeshX*MeshY))));
            SimuNodes.at(i).setNbhd(WHICHDIR::DirD, (i+MeshX*MeshY >= Num_Nodes ? &(SimuNodes.at(i)) : &(SimuNodes.at(i+MeshX*MeshY))));

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
                //down
                if (i>=0 && i<MeshX*MeshY){
                    SimuNodes.at(i).setNbhd(WHICHDIR::DirU, &(SimuNodes.at(i+MeshX*MeshY*(MeshZ-1))));
                }
                //up
                if (i>=MeshX*MeshY*(MeshZ-1) && i<MeshX*MeshY*MeshZ){
                    SimuNodes.at(i).setNbhd(WHICHDIR::DirD, &(SimuNodes.at(i-MeshX*MeshY*(MeshZ-1))));
                }
            }

            if (whichBOUNDCOND == BOUNDCOND::BndCONST){
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
            
            if(whichBOUNDCOND == BOUNDCOND::BndADIABATIC){}
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
        case WHICHPARA::CUSTOM:
            for(auto &node : SimuNodes)node.Cust_Node.deletEntry(_Index);
            break;                    
        default:
            break;
        }
    }
    
    /*************************************************************/

    // get value of whichpara at Index's entry
    std::vector<double> getMeshProp(WHICHPARA whichpara, int Index){
        std::vector<double> result(Num_Nodes, 0.0);
        if (getNum_Ent(whichpara) != 0 &&Index<getNum_Ent(whichpara))
        #pragma omp parallel for
            for (int i = 0; i < Num_Nodes;i++){
                result.at(i) = (SimuNodes.at(i).getVal(whichpara,Index));
            }
        else throw std::out_of_range("Index out of range");
        return result;
    }

    std::vector<double> getUni_Prop(WHICHPARA whichpara){
        std::vector<double> result;
        for (auto node : SimuNodes){
            double sum = 0;
            for (auto val : node.getVal(whichpara)){
                sum += val*val;
            }
            sum = sqrt(sum);
            threshold(sum, 0.0001, 0.9999);
            result.push_back(sum);
        }
        return result;
    }

    /*************************************************************/

    void threshold(double &val, double min, double max){
        val>max ? val = max : (val<min ? val = min : val = val);
    }


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

    void updateNodeCon(int where, int Index, double _phs){
        SimuNodes.at(where).Con_Node.updateVal(Index, _phs);
    }

    void updateNodeCon(std::vector<int> where, int Index, double _phs){
        updateNodeCon(transCoord(where), Index, _phs);
    }

    /*************************************************************/

    void updateNodePhs(int where, int Index, double _phs){
        SimuNodes.at(where).Phs_Node.updateVal(Index, _phs);
    }

    void updateNodePhs(std::vector<int> where, int Index, double _phs){
        updateNodePhs(transCoord(where), Index, _phs);
    }


    /*************************************************************/

    /**/
    void Laplacian(STENCILE whichSTNCL, WHICHPARA whichpara){
        double result = 0;
        double dx = StepX;
        double dy = StepY;
        double dz = StepZ;
        if (whichSTNCL == STENCILE::StencFIVE){
        #pragma omp parallel for collapse(2)
        for (auto &node : SimuNodes){
            for (int Index = 0; Index < SimuNodes.at(0).getNum_Ent(whichpara); ++Index){
                    double c = node.getVal(whichpara,Index);
                    double f = node.getNbhd(WHICHDIR::DirF)->getVal(whichpara,Index);
                    double b = node.getNbhd(WHICHDIR::DirB)->getVal(whichpara,Index);
                    double l = node.getNbhd(WHICHDIR::DirL)->getVal(whichpara,Index);
                    double r = node.getNbhd(WHICHDIR::DirR)->getVal(whichpara,Index);

                    result = ((f+b+l+r-4*c)/(dx*dy*dz));
                    if (whichpara == WHICHPARA::PHSFRAC)
                        node.Phs_Node.updateLap(Index, result);
                    if (whichpara == WHICHPARA::CON)
                        node.Con_Node.updateLap(Index, result);
                    if (whichpara == WHICHPARA::CUSTOM)
                        node.Cust_Node.updateLap(Index, result);
                }
            }
        }
        return;
    }

    void Laplacian(WHICHPARA whichpara){
        Laplacian(STENCILE::StencFIVE, whichpara);
        return;
    }
    /*************************************************************/
    void showGlobalInfo(); 
    void showNodesProp(WHICHPARA which, int Index); 

    void outVTKFilehead(std::string _dirname, int istep);
    void outVTKAve(std::string _dirname, WHICHPARA whichpara, int istep);
    void outVTKAll(std::string _dirname, WHICHPARA whichpara, int istep);

    void outVTK(std::string _dirname, int istep){
        outVTKFilehead(_dirname, istep);
        outVTKAve(_dirname, WHICHPARA::CON, istep);
        outVTKAve(_dirname, WHICHPARA::PHSFRAC, istep);
    }
};

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
void SimulationMesh::showNodesProp(WHICHPARA which, int Index){ //which para, Index of para
    if (which == WHICHPARA::CON)std::cout<<"Concentration of "<<SimuNodes.at(0).Con_Node.Entrys.at(Index).Element<<"\n";
    if (which == WHICHPARA::PHSFRAC)std::cout<<"Order Parameter of "<<Index<<" grain\n";
    if (which == WHICHPARA::CUSTOM)std::cout<<"Custom value \n";
    for (long i = 0; i < Num_Nodes; i++){
        std::cout<<std::fixed<<std::setprecision(10)<<
            SimuNodes.at(i).getVal(which,Index)
            // SimuNodes.at(i).getLap(WHICHPARA::CON, Index)
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


    char *filename = new char[_dirname.length()+1];
    std::strcpy(filename, _dirname.c_str());

    std::sprintf(filename, "%s/time_%04d.vtk", filename, istep);

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
    for (int i = 0;i<MeshX;i++)
        for (int j = 0;j<MeshY;j++)
            for (int k = 0; k < MeshZ; k++){
                outfile<<i*StepLength.at(0)<<"   "<<j*StepLength.at(1)<<"   "<<k*StepLength.at(2)<<"\n";
            }
    outfile<<"POINT_DATA "<<Num_Nodes<<"\n";

    outfile.close();
}

inline void SimulationMesh::outVTKAve(std::string _dirname, WHICHPARA whichpara, int istep){

    char *filename = new char[_dirname.length()+1];
    std::strcpy(filename, _dirname.c_str());
    std::sprintf(filename, "%s/time_%04d.vtk", filename, istep);

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
    default:
        break;
    }
    std::vector<double> normed_phs(getUni_Prop(whichpara));
    outfile<<"SCALARS "<<varname<<"  float  1\n";
    outfile<<"LOOKUP_TABLE default\n";
    for (int i = 0;i<Num_Nodes;i++)
        outfile<<normed_phs.at(i)<<"\n";

    outfile.close();

}

inline void SimulationMesh::outVTKAll(std::string _dirname, WHICHPARA whichpara, int istep){

    char *filename = new char[_dirname.length()+1];
    std::strcpy(filename, _dirname.c_str());
    std::sprintf(filename, "%s/time_%04d.vtk", filename, istep);

    std::ofstream outfile;
    outfile.open(filename, std::ios::app);

    char varname[64];

    std::vector<std::vector<double>> meshphs;
    meshphs.reserve(getNum_Ent(whichpara));
    for (int num = 0; num < getNum_Ent(whichpara); num++){
        meshphs.push_back(getMeshProp(whichpara, num));
    }

    for (int num = 0; num < getNum_Ent(whichpara); num++){
        switch (whichpara){
        case WHICHPARA::PHSFRAC:
            std::sprintf(varname, "PHSFRAC_%01d", num);
            break;
        case WHICHPARA::CON:
            std::sprintf(varname, "CON__%01d", num);
            break;
        case WHICHPARA::CUSTOM:
            std::sprintf(varname, "CUST__%01d", num);
            break;
        default:
            break;
        }
        outfile<<"SCALARS "<<varname<<"  float  1\n";
        outfile<<"LOOKUP_TABLE default\n";
        for (int i = 0;i<Num_Nodes;i++)
            outfile<<meshphs.at(num).at(i)<<"\n";
    }
    outfile.close();
}


#endif