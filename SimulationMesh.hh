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
    private:
        std::vector<double> Dimension{64, 64, 1};
        double& BoxX = Dimension.at(0);
        double& BoxY = Dimension.at(1);
        double& BoxZ = Dimension.at(2);

        std::vector<double> StepLength{1,1,1};
        // double& dx = StepLength.at(0);
        // double& dy = StepLength.at(1);
        // double& dz = StepLength.at(2);

        std::vector<MeshNode> SimuNodes;
    public:

        SimulationMesh(){} // initial with default size, without SimuNodes.at(i)s

        SimulationMesh(std::vector<double> SizeInfo, MeshNode Node){ // initial with size and SimuNodes.at(i)s
            Dimension = SizeInfo;
            fillNodes(Node);
            bindBoundary(BOUNDCOND::PERIODIC);
        }
        SimulationMesh(MeshNode Node){ // initial with SimuNodes.at(i)s and default size
            fillNodes(Node);
            bindBoundary(BOUNDCOND::PERIODIC);
        }
        ~SimulationMesh(){};

/*************************************************************/

        void fillNodes(MeshNode Nodes){ // fill mesh with SimuNodes.at(i)s 
            for(double i = 0; i < Dimension.at(0)*Dimension.at(1)*Dimension.at(2); i++){
                SimuNodes.push_back(Nodes);
            }
        }

        void bindBoundary(double whichBOUNDCOND){
            for(double i = 0; i < getNum_Nodes(); i++ ){
                SimuNodes.at(i).whichdir.Forward = (i-BoxY < 0? &(SimuNodes.at(i)): &(SimuNodes.at(i-BoxY)) );
                SimuNodes.at(i).whichdir.Backward = (i+BoxY >= getNum_Nodes()? &(SimuNodes.at(i)) : &(SimuNodes.at(i+BoxY)));
                SimuNodes.at(i).whichdir.Left =(i-1 < 0? &(SimuNodes.at(i)) : &(SimuNodes.at(i-1)));
                SimuNodes.at(i).whichdir.Right = (i+1 >= getNum_Nodes()? &(SimuNodes.at(i)) : &(SimuNodes.at(i+1)));
                SimuNodes.at(i).whichdir.Down= (i-BoxX*BoxY < 0? &(SimuNodes.at(i)): &(SimuNodes.at(i-BoxX*BoxY)));
                SimuNodes.at(i).whichdir.Up= (i+BoxX*BoxY >= getNum_Nodes()? &(SimuNodes.at(i)): &(SimuNodes.at(i+BoxX*BoxY)));

                if(whichBOUNDCOND == BOUNDCOND::PERIODIC){
                    if((int)i%(int)BoxX == 0){
                        SimuNodes.at(i).whichdir.Left = &(SimuNodes.at(i+(BoxX-1)));
                    }
                    if((int)i%(int)BoxX == (BoxX-1)){
                        SimuNodes.at(i).whichdir.Right = &(SimuNodes.at(i-(BoxX-1)));
                    }
                    if(((int)i%(int)(BoxX*BoxY))>=0 && ((int)i%(int)(BoxX*BoxY))<BoxX){
                        SimuNodes.at(i).whichdir.Forward = &(SimuNodes.at(i+BoxX*(BoxY-1)));
                    }
                    if(((int)i%(int)(BoxX*BoxY))>=BoxX*(BoxY-1) && ((int)i%(int)(BoxX*BoxY))<BoxX*BoxY){
                        SimuNodes.at(i).whichdir.Backward = &(SimuNodes.at(i-BoxX*(BoxY-1)));
                    }
                    if(i>=0 && i<BoxX*BoxY){
                        SimuNodes.at(i).whichdir.Down = &(SimuNodes.at(i+BoxX*BoxY*(BoxZ-1)));
                    }
                    if(i>=BoxX*BoxY*(BoxZ-1) && i<BoxX*BoxY*BoxZ){
                        SimuNodes.at(i).whichdir.Up = &(SimuNodes.at(i-BoxX*BoxY*(BoxZ-1)));
                    }
                }
//              /**/
            }
        }

/*************************************************************/

        std::vector<double> getMeshProp(unsigned which, unsigned index){
            std::vector<double> result(getNum_Nodes(),0.0);
            for(unsigned i = 0; i < getNum_Nodes();i++){
                result.at(i) = (SimuNodes.at(i).getProp(which).at(index));
            }
            return result;
        }

        std::vector<double> getMeshCon(ELEMENT _element){
            std::vector<double> meshcon(getNum_Nodes(),0.0);
            for(unsigned i = 0; i < getNum_Nodes();i++){
                meshcon.at(i) = (SimuNodes.at(i).Con_Node.getCon(_element));
            }
            return meshcon;
        }

        double getDim(const unsigned which){
            return Dimension.at(which);
        }

        double getNum_Nodes(){ // return the number of SimuNodes.at(i)s in mesh
            return SimuNodes.size();
        }

        double getStepLength(const int which){
            return StepLength.at(which);
        }

        unsigned getNum_Prop(WHICHPARA which){
            return SimuNodes.at(0).getNum(which);
        }

        std::vector<double> getUni_Prop(WHICHPARA whichpara){
            std::vector<double> result;
            for(auto node : SimuNodes){
                double sum = 0;
                for(auto val : node.getProp(whichpara)){
                    sum += val*val;
                }
                sum = sqrt(sum) / node.getNum(whichpara);
                result.push_back(sum);
            }
            return result;
        }


        std::vector<double> transCoord(double where){
            if(where<getNum_Nodes()){       
                std::vector<double>coord;
                coord.at(0) = (int)where%(int)BoxX;
                coord.at(1) = (int)((where-coord.at(0)) / BoxX)%(int)BoxY;
                coord.at(2) = (int)(where-(coord.at(0)+coord.at(1)*BoxX)/(BoxX*BoxY));
                return coord;
            }
            throw std::out_of_range("Not in mesh");
            return {};
        }

        double transCoord(std::vector<double> where){
            if(where.at(0)<BoxX&&where.at(1)<BoxY&&where.at(2)<BoxZ)
            return where.at(0)+where.at(1)*BoxX+where.at(2)*BoxX*BoxY;
            throw std::out_of_range("Not in mesh");
            return 0;
        }
        
        // find SimuNodes.at(i) in the mesh according to the coordinates
        MeshNode& findNode(double X, double Y, double Z){ 
            if(X<Dimension.at(0)&&Y<Dimension.at(1)&&Z<Dimension.at(2) && !(X<0) &&!(Y<0) &&!(Z<0)){
                return SimuNodes.at(X+Y*Dimension.at(0)+Z*Dimension.at(0)*Dimension.at(1));
            }
            else throw std::invalid_argument("Index Not in Mesh");
        }

        MeshNode& findNode(double AbsCoord){
            if(AbsCoord<getNum_Nodes()){
                return SimuNodes.at(AbsCoord);
            }
            else throw std::invalid_argument("Index Not in Mesh");
        }

        void updateCustomValue(std::vector<double> CustomValue){
            for(int i = 0; i < getNum_Nodes(); i++){
                SimuNodes.at(i).Custom_Value = CustomValue.at(i);
            }
        }

/*************************************************************/

        void updateMeshCon(std::vector<double> _val){
            for(int i = 0; i < getNum_Nodes(); i++){
                if(SimuNodes.at(i).getNum(WHICHPARA::CON) == 1)
                SimuNodes.at(i).Con_Node.updateEntry(_val.at(i));
                else throw std::invalid_argument("Exist more than one element");
            }
        }

        void updateMeshCon(ELEMENT _element,std::vector<double> _val){
            for(int i = 0; i < getNum_Nodes(); i++){
                SimuNodes.at(i).Con_Node.updateEntry(_element,_val.at(i));
            }
        }

        void updateNodeCon(unsigned where, double _con){
            if(SimuNodes.at(where).getNum(WHICHPARA::CON) == 1)
            SimuNodes.at(where).Con_Node.updateEntry(_con);
            else throw std::invalid_argument("Exist more than one element");
        }

        void updateNodeCon(unsigned where, ELEMENT _element, double _con){
            SimuNodes.at(where).Con_Node.updateEntry(_element,_con);
        } 

        void updateNodeCon(std::vector<double> where, double _con){
            updateNodeCon(transCoord(where) , _con);
        }

        void updateNodeCon(std::vector<double> where, ELEMENT _element, double _con){
            updateNodeCon(transCoord(where),_element,_con);
        }

/*************************************************************/
        void updateMeshPhs(unsigned index, std::vector<double> _val){
            for(unsigned i = 0; i < getNum_Nodes(); i++){
                SimuNodes.at(i).Phs_Node.updateEntry(index,_val.at(i));
            }
        }

        void updateNodePhs(unsigned where,unsigned index, double _phs){
            SimuNodes.at(where).Phs_Node.updateEntry(index,_phs);
        }

        void updateNodePhs(std::vector<double> where, unsigned index, double _phs){
            updateNodePhs(transCoord(where),index,_phs);
        }



/*************************************************************/
        void showGlobalInfo(); // show the basic information of the mesh
        void showNodesProp(unsigned which, unsigned index); // show one of the properties of the SimuNodes.at(i)s inside the mesh
        void write_vtk_grid_values(int istep,std::vector<double> data);
        void outFile(int istep, std::vector<double> data);

/*************************************************************/

        std::vector<double> Laplacian (int whichSTNCL, int whichPara){
            std::vector<double> result(getNum_Nodes(),0.0);
            double dx = getStepLength(WHICHDIM::X);
            double dy = getStepLength(WHICHDIM::Y);
            double dz = getStepLength(WHICHDIM::Z);
            if(whichSTNCL = STENCILE::FIVEPOINT){
            #pragma omp parallel for
                for(int i = 0; i < getNum_Nodes(); i++){
                    
                    double c = findNode(i).getProp(whichPara).at(0);
                    double f = findNode(i).whichdir.Forward->getProp(whichPara).at(0);
                    double b = findNode(i).whichdir.Backward->getProp(whichPara).at(0);
                    double l = findNode(i).whichdir.Left->getProp(whichPara).at(0);
                    double r = findNode(i).whichdir.Right->getProp(whichPara).at(0);
    
                    result.at(i) = ((f+b+l+r-4*c)/(dx*dy*dz));
                }
            }
            return result;
        }
};

void SimulationMesh::showGlobalInfo(){
    std::cout<<"SimulationMesh Properties:"<<std::endl;
    std::cout<<"Mesh Size:\t\t"<<BoxX<<"\u0078"<<BoxY<<"\u0078"<<BoxZ<<std::endl;
    std::cout<<"Number of Nodes:\t"<<getNum_Nodes()<<std::endl;
    std::cout<<"Elements:\t\t";
    for(auto elements : SimuNodes.at(0).Con_Node.getElementList())
    std::cout<<elements<<" ";
    std::cout<<"\nNumber of Phase:\t"<<SimuNodes.at(0).getNum(WHICHPARA::PHSFRAC);

    std::cout<<"\n-----------------------------------------------------------------\n"<<std::endl;
}

void SimulationMesh::showNodesProp(unsigned which, unsigned index){ //which para, index of para
    if(which == WHICHPARA::CON)std::cout<<"Concentration of "<<SimuNodes.at(0).Con_Node.getElementList().at(index)<<"\n";
    if(which == WHICHPARA::PHSFRAC)std::cout<<"Order Parameter of "<<index<<" grain\n";
    if(which == WHICHPARA::CUSTOM)std::cout<<"Custom value \n";
    for(long i = 0; i < getNum_Nodes(); i++){
        std::cout<<std::fixed<<std::setprecision(10)<<SimuNodes.at(i).getProp(which).at(index)<<" ";
        if(i%(int)BoxX==BoxX-1)std::cout<<"\n";
        if(i%(int)(BoxX*BoxY) == BoxX*BoxY-1)std::cout<<"\n";
    }
    std::cout<<"-----------------------------------------------------------------\n"<<std::endl;
}

void SimulationMesh::write_vtk_grid_values(int istep,std::vector<double> data)
{

	char filename[128];
	sprintf(filename,"D:\\Developing\\C++\\PhaseFieldModelling\\CaseStudy_3\\output\\Result\\time_%04d.vtk", istep);

	std::ofstream outfile;
	outfile.open(filename);
	outfile<<"# vtk DataFile Version 2.0\n";
	outfile<<"time_10.vtk\n";
	outfile<<"ASCII\n";
	outfile<<"DATASET STRUCTURED_GRID\n";

	outfile<<"DIMENSIONS "<<BoxX<<"  "<<BoxY<<"  "<<BoxZ<<"\n";
	outfile<<"POINTS "<<getNum_Nodes()<<"   float\n";
	double dumx,dumy,dumz;
	for(int i=0;i<BoxX;i++)
		for(int j=0;j<BoxY;j++)
		{
			dumx = i*StepLength.at(0); dumy = j*StepLength.at(1); dumz = 0.0;

			outfile<<dumx<<"   "<<dumy<<"   "<<dumz<<"\n";
		}

	outfile<<"POINT_DATA "<<getNum_Nodes()<<"\n";
	outfile<<"SCALARS PHSFRAC  float  1\n";
	outfile<<"LOOKUP_TABLE default\n";

	for(int i=0;i<BoxX*BoxY;i++)
		outfile<<data.at(i)<<"\n";

	outfile.close();
	
}

void SimulationMesh::outFile(int istep, std::vector<double> data){
    write_vtk_grid_values(istep,data);
}


#endif