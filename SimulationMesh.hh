#pragma once
#ifndef PHASE_SIMULATION_MESH
#define PHASE_SIMULATION_MESH

#include <iostream>
#include <vector>
#include <fstream>
#include "MeshNode.hh"

enum WHICHDIM{X,Y,Z};
enum BOUNDCOND{PERIODIC,CONST,ADIABATIC};
enum STENCILE{FIVEPOINT=5,NINEPOINT=9};


class SimulationMesh{
    public:
        std::vector<int> Dimension{64, 64, 1};
        int& MeshX = Dimension.at(0);
        int& MeshY = Dimension.at(1);
        int& MeshZ = Dimension.at(2);

        std::vector<double> StepLength{1,1,1};
        double& StepX = StepLength.at(0);
        double& StepY = StepLength.at(1);
        double& StepZ = StepLength.at(2);

        int Num_Nodes;

        std::vector<MeshNode> SimuNodes;

/*************************************************************/

        SimulationMesh(){} // initial with default size, without (*this)(i)s

        SimulationMesh(std::vector<int> SizeInfo,std::vector<double> StepInfo, MeshNode Node){ // initial with size and (*this)(i)s
            Dimension = SizeInfo;
            StepLength = StepInfo;
            fillNodes(Node);
            Num_Nodes = SimuNodes.size();
            bindBoundary(BOUNDCOND::PERIODIC);
        }

        SimulationMesh(std::vector<double> StepInfo,MeshNode Node):SimulationMesh({64,64,1},StepInfo,Node){}
        SimulationMesh(MeshNode Node):SimulationMesh({1,1,1},Node){} // initial with (*this)(i)s and default size

        ~SimulationMesh(){
            Dimension.clear();
            StepLength.clear();
            SimuNodes.clear();
        };

        MeshNode& operator()(const int where){
            if(where<Num_Nodes){
                return SimuNodes.at(where);
            }
            else throw std::invalid_argument("Index Not in Mesh");
        }

        MeshNode& operator()(int X, int Y, int Z){ 
            if(X<Dimension.at(0)&&Y<Dimension.at(1)&&Z<Dimension.at(2) && !(X<0) &&!(Y<0) &&!(Z<0)){
                return SimuNodes.at(X+Y*Dimension.at(0)+Z*Dimension.at(0)*Dimension.at(1));
            }
            else throw std::invalid_argument("Index Not in Mesh");
        }


        void addEntry(WHICHPARA whichpara, int num){
            switch (whichpara)
            {
            case WHICHPARA::CON :
                for(auto & node : SimuNodes)node.Con_Node.addEntry(num);
                break;
            case WHICHPARA::PHSFRAC :
                for(auto & node : SimuNodes)node.Phs_Node.addEntry(num);
                break;
            default:
                break;
            }
        }

        /*************************************************************/

        void fillNodes(MeshNode Nodes){ // fill mesh with (*this)(i)s 
            SimuNodes.reserve(MeshX*MeshY*MeshZ);
            for(int i = 0; i < MeshX*MeshY*MeshZ; i++){
                SimuNodes.push_back(Nodes);
            }
        }

        void bindBoundary(BOUNDCOND whichBOUNDCOND){
            for(int i = 0; i < Num_Nodes; i++ ){
                (*this)(i).Forward = (i-MeshY < 0? &((*this)(i)): &((*this)(i-MeshY)) );
                (*this)(i).Backward = (i+MeshY >= Num_Nodes? &((*this)(i)) : &((*this)(i+MeshY)));
                (*this)(i).Left =(i-1 < 0? &((*this)(i)) : &((*this)(i-1)));
                (*this)(i).Right = (i+1 >= Num_Nodes? &((*this)(i)) : &((*this)(i+1)));
                (*this)(i).Down= (i-MeshX*MeshY < 0? &((*this)(i)): &((*this)(i-MeshX*MeshY)));
                (*this)(i).Up= (i+MeshX*MeshY >= Num_Nodes? &((*this)(i)): &((*this)(i+MeshX*MeshY)));

                if(whichBOUNDCOND == BOUNDCOND::PERIODIC){
                    if(i%MeshX == 0){
                        (*this)(i).Left = &((*this)(i+(MeshX-1)));
                    }
                    if(i%MeshX == (MeshX-1)){
                        (*this)(i).Right = &((*this)(i-(MeshX-1)));
                    }
                    if((i%(MeshX*MeshY))>=0 && (i%(MeshX*MeshY))<MeshX){
                        (*this)(i).Forward = &((*this)(i+MeshX*(MeshY-1)));
                    }
                    if((i%(MeshX*MeshY))>=MeshX*(MeshY-1) && (i%(MeshX*MeshY))<MeshX*MeshY){
                        (*this)(i).Backward = &((*this)(i-MeshX*(MeshY-1)));
                    }
                    if(i>=0 && i<MeshX*MeshY){
                        (*this)(i).Down = &((*this)(i+MeshX*MeshY*(MeshZ-1)));
                    }
                    if(i>=MeshX*MeshY*(MeshZ-1) && i<MeshX*MeshY*MeshZ){
                        (*this)(i).Up = &((*this)(i-MeshX*MeshY*(MeshZ-1)));
                    }
                }
//              /**/
            }
        }

        /*************************************************************/

        std::vector<double> getMeshProp(WHICHPARA whichpara, int index){
            std::vector<double> result(Num_Nodes,0.0);
            #pragma parallel for
            if(index<getNum_Prop(whichpara))
            for(int i = 0; i < Num_Nodes;i++){
                result.at(i) = ((*this)(i).getProp(whichpara).at(index));
            }
            else throw std::out_of_range("index out of range");
            return result;
        }

        std::vector<double> getUni_Prop(WHICHPARA whichpara){
            std::vector<double> result;
            for(auto node : SimuNodes){
                double sum = 0;
                for(auto val : node.getProp(whichpara)){
                    sum += val*val;
                }
                sum = sqrt(sum);
                threshold(sum,0.0001,0.9999);
                result.push_back(sum);
            }
            return result;
        }

        void threshold(double &val, double min, double max){
            val>max?val = max:(val<min?val = min:val = val);
        }

        /*************************************************************/

        std::vector<int> transCoord(int where){
            if(where<Num_Nodes){       
                std::vector<int>coord;
                coord.at(0) = where%MeshX;
                coord.at(1) = ((where-coord.at(0)) / MeshX)%MeshY;
                coord.at(2) = (where-(coord.at(0)+coord.at(1)*MeshX)/(MeshX*MeshY));
                return coord;
            }
            throw std::out_of_range("Not in mesh");
            return {};
        }

        int transCoord(std::vector<int> where){
            if(where.at(0)<MeshX&&where.at(1)<MeshY&&where.at(2)<MeshZ)
            return where.at(0)+where.at(1)*MeshX+where.at(2)*MeshX*MeshY;
            throw std::out_of_range("Not in mesh");
            return 0;
        }

        int getNum_Prop(WHICHPARA which){
            return (*this)(0).getNum_Ent(which);
        }
        
        void updateCustValue(int _where,int _index, double _val){
            (*this)(_where).Cust_Node.CustVal.at(_index) = _val;
        }

        /*************************************************************/

        void updateNodeCon(int where, double _con){
            if((*this)(where).getNum_Ent(WHICHPARA::CON) == 1)
            (*this)(where).Con_Node.updateEntry(0,_con);
            else throw std::invalid_argument("Exist more than one element");
        }
        
        void updateNodeCon(std::vector<int> where, double _con){
            updateNodeCon(transCoord(where) , _con);
        }

        /*******/

        void updateNodeCon(int where, ELEMENT _element, double _con){
            (*this)(where).Con_Node.updateEntry(_element,_con);
        } 

        void updateNodeCon(std::vector<int> where, ELEMENT _element, double _con){
            updateNodeCon(transCoord(where),_element,_con);
        }

        /*************************************************************/
        void updateMeshPhs(int index, std::vector<double> _val){
            for(int i = 0; i < Num_Nodes; i++){
                (*this)(i).Phs_Node.updateEntry(index,_val.at(i));
            }
        }

        void updateNodePhs(int where,int index, double _phs){
            (*this)(where).Phs_Node.updateEntry(index,_phs);
        }

        void updateNodePhs(std::vector<int> where, int index, double _phs){
            updateNodePhs(transCoord(where),index,_phs);
        }

/*************************************************************/
        void showGlobalInfo(); // show the basic information of the mesh
        void showNodesProp(WHICHPARA which, int index); // show one of the properties of the (*this)(i)s inside the mesh
        void outFile(int istep);

/*************************************************************/

/**/    void Laplacian (STENCILE whichSTNCL, WHICHPARA whichpara){
            double result = 0;
            double dx = StepX;
            double dy = StepY;
            double dz = StepZ;
            if(whichSTNCL = STENCILE::FIVEPOINT){
                #pragma omp parallel for
                for(int index = 0; index < (*this)(0).getNum_Ent(whichpara); ++index){
                    for(auto &node : SimuNodes){
                        double c = node.getProp(whichpara).at(index);
                        double f = node.Forward->getProp(whichpara).at(index);
                        double b = node.Backward->getProp(whichpara).at(index);
                        double l = node.Left->getProp(whichpara).at(index);
                        double r = node.Right->getProp(whichpara).at(index);

                        result = ((f+b+l+r-4*c)/(dx*dy*dz));
                        if(whichpara == WHICHPARA::PHSFRAC)
                            node.Phs_Node.updateLap(index,result);
                        if(whichpara == WHICHPARA::CON)
                            node.Con_Node.updateLap(index,result);
                        if(whichpara == WHICHPARA::CUSTOM)
                            node.Cust_Node.updateLap(index,result);
                    }
                }
            }
            return;
        }

        void Laplacian (WHICHPARA whichpara){
            Laplacian(STENCILE::FIVEPOINT,whichpara);
            return ;
        }
};

void SimulationMesh::showGlobalInfo(){
    std::cout<<"SimulationMesh Properties:\n";
    std::cout<<"Mesh Size:\t\t"<<MeshX<<"\u0078"<<MeshY<<"\u0078"<<MeshZ<<"\n";
    std::cout<<"Elements:\t\t";
    for(auto ent : (*this)(0).Con_Node.Entrys)
    std::cout<<ent.Element<<" ";
    std::cout<<"\nNumber of Phase:\t"<<(*this)(0).getNum_Ent(WHICHPARA::PHSFRAC);

    std::cout<<"\n-----------------------------------------------------------------\n"<<std::endl;
}

void SimulationMesh::showNodesProp(WHICHPARA which, int index){ //which para, index of para
    if(which == WHICHPARA::CON)std::cout<<"Concentration of "<<(*this)(0).Con_Node.Entrys.at(index).Element<<"\n";
    if(which == WHICHPARA::PHSFRAC)std::cout<<"Order Parameter of "<<index<<" grain\n";
    if(which == WHICHPARA::CUSTOM)std::cout<<"Custom value \n";
    for(long i = 0; i < Num_Nodes; i++){
        std::cout<<std::fixed<<std::setprecision(10)<<
            // (*this)(i).getProp(which).at(index)
            (*this)(i).getLap(WHICHPARA::CON,index)
        <<" ";
        if(i%MeshX==MeshX-1)std::cout<<"\n";
        if(i%(MeshX*MeshY) == MeshX*MeshY-1)std::cout<<"\n";
    }
    std::cout<<"-----------------------------------------------------------------\n"<<std::endl;
}

inline void SimulationMesh::outFile(int istep){

	char filename[128];
	sprintf(filename,"../output/Result_1/time_%04d.vtk", istep);

	std::ofstream outfile;
	outfile.open(filename);
	outfile<<"# vtk DataFile Version 2.0\n";
	outfile<<"time_10.vtk\n";
	outfile<<"ASCII\n";
	outfile<<"DATASET STRUCTURED_GRID\n";

	outfile<<"DIMENSIONS "<<MeshX<<"  "<<MeshY<<"  "<<MeshZ<<"\n";
	outfile<<"POINTS "<<Num_Nodes<<"   float\n";
	double dumx,dumy,dumz;
	for(int i=0;i<MeshX;i++)
		for(int j=0;j<MeshY;j++)
            for(int k = 0; k < MeshZ; k++)
		    {
		    	dumx = i*StepLength.at(0); dumy = j*StepLength.at(1); dumz = k*StepLength.at(2);
    
		    	outfile<<dumx<<"   "<<dumy<<"   "<<dumz<<"\n";
		    }

	outfile<<"POINT_DATA "<<Num_Nodes<<"\n";

    char varname[64];
    // int phs_num = getNum_Prop(WHICHPARA::PHSFRAC);
    // int con_num = getNum_Prop(WHICHPARA::CON);


        std::vector<std::vector<double>> meshcon;

    for(int num = 0; num < getNum_Prop(WHICHPARA::CON); num ++){
        meshcon.push_back(getMeshProp(WHICHPARA::CON,num));
    }

	for(int num = 0; num < getNum_Prop(WHICHPARA::CON); num++){
        sprintf(varname,"CON_%01d",num);
	    outfile<<"SCALARS "<<varname<<"  float  1\n";
	    outfile<<"LOOKUP_TABLE default\n";

	    for(int i=0;i<Num_Nodes;i++)
	    	outfile<<meshcon.at(num).at(i)<<"\n";
    }

    std::vector<std::vector<double>> meshphs;

    meshphs.reserve(getNum_Prop(WHICHPARA::PHSFRAC));
    for(int num = 0; num < getNum_Prop(WHICHPARA::PHSFRAC); num ++){
        meshphs.push_back(getMeshProp(WHICHPARA::PHSFRAC,num));
    }

    // # parallel for
    for(int num = 0; num < getNum_Prop(WHICHPARA::PHSFRAC); num++){
        sprintf(varname,"PHSFRAC_%01d",num);
	    outfile<<"SCALARS "<<varname<<"  float  1\n";
	    outfile<<"LOOKUP_TABLE default\n";
	    for(int i=0;i<Num_Nodes;i++)
	    	outfile<<meshphs.at(num).at(i)<<"\n";
    }

        std::vector<double> normed_phs(getUni_Prop(WHICHPARA::PHSFRAC));
        sprintf(varname,"PHSFRAC_NORMED");
	    outfile<<"SCALARS "<<varname<<"  float  1\n";
	    outfile<<"LOOKUP_TABLE default\n";
	    for(int i=0;i<Num_Nodes;i++)
	    	outfile<<normed_phs.at(i)<<"\n";
    
        std::vector<double> cust = getMeshProp(WHICHPARA::CUSTOM,0);
        sprintf(varname,"CUSTOM");
	    outfile<<"SCALARS "<<varname<<"  float  1\n";
	    outfile<<"LOOKUP_TABLE default\n";
	    for(int i=0;i<Num_Nodes;i++)
	    	outfile<<cust.at(i)<<"\n";

    meshphs.clear();
    meshcon.clear();

    outfile.close();

}

#endif