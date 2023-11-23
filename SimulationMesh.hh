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
        std::vector<int> Dimension{64, 64, 1};
        int& BoxX = Dimension.at(0);
        int& BoxY = Dimension.at(1);
        int& BoxZ = Dimension.at(2);

        std::vector<double> StepLength{1,1,1};

        std::vector<MeshNode> SimuNodes;
        int Num_Nodes;
    public:

        SimulationMesh(){} // initial with default size, without SimuNodes.at(i)s

        SimulationMesh(std::vector<int> SizeInfo,std::vector<double> StepInfo, MeshNode Node){ // initial with size and SimuNodes.at(i)s
            Dimension = SizeInfo;
            StepLength = StepInfo;
            fillNodes(Node);
            Num_Nodes = SimuNodes.size();
            bindBoundary(BOUNDCOND::PERIODIC);
        }

        SimulationMesh(std::vector<double> StepInfo,MeshNode Node):SimulationMesh({64,64,1},StepInfo,Node){}
        SimulationMesh(MeshNode Node):SimulationMesh({1,1,1},Node){} // initial with SimuNodes.at(i)s and default size

        ~SimulationMesh(){
            Dimension.clear();
            StepLength.clear();
            SimuNodes.clear();
        };

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

        void fillNodes(MeshNode Nodes){ // fill mesh with SimuNodes.at(i)s 
            SimuNodes.reserve(BoxX*BoxY*BoxZ);
            for(int i = 0; i < BoxX*BoxY*BoxZ; i++){
                SimuNodes.push_back(Nodes);
            }
        }

        void bindBoundary(BOUNDCOND whichBOUNDCOND){
            for(int i = 0; i < Num_Nodes; i++ ){
                SimuNodes.at(i).Forward = (i-BoxY < 0? &(SimuNodes.at(i)): &(SimuNodes.at(i-BoxY)) );
                SimuNodes.at(i).Backward = (i+BoxY >= Num_Nodes? &(SimuNodes.at(i)) : &(SimuNodes.at(i+BoxY)));
                SimuNodes.at(i).Left =(i-1 < 0? &(SimuNodes.at(i)) : &(SimuNodes.at(i-1)));
                SimuNodes.at(i).Right = (i+1 >= Num_Nodes? &(SimuNodes.at(i)) : &(SimuNodes.at(i+1)));
                SimuNodes.at(i).Down= (i-BoxX*BoxY < 0? &(SimuNodes.at(i)): &(SimuNodes.at(i-BoxX*BoxY)));
                SimuNodes.at(i).Up= (i+BoxX*BoxY >= Num_Nodes? &(SimuNodes.at(i)): &(SimuNodes.at(i+BoxX*BoxY)));

                if(whichBOUNDCOND == BOUNDCOND::PERIODIC){
                    if(i%BoxX == 0){
                        SimuNodes.at(i).Left = &(SimuNodes.at(i+(BoxX-1)));
                    }
                    if(i%BoxX == (BoxX-1)){
                        SimuNodes.at(i).Right = &(SimuNodes.at(i-(BoxX-1)));
                    }
                    if((i%(BoxX*BoxY))>=0 && (i%(BoxX*BoxY))<BoxX){
                        SimuNodes.at(i).Forward = &(SimuNodes.at(i+BoxX*(BoxY-1)));
                    }
                    if((i%(BoxX*BoxY))>=BoxX*(BoxY-1) && (i%(BoxX*BoxY))<BoxX*BoxY){
                        SimuNodes.at(i).Backward = &(SimuNodes.at(i-BoxX*(BoxY-1)));
                    }
                    if(i>=0 && i<BoxX*BoxY){
                        SimuNodes.at(i).Down = &(SimuNodes.at(i+BoxX*BoxY*(BoxZ-1)));
                    }
                    if(i>=BoxX*BoxY*(BoxZ-1) && i<BoxX*BoxY*BoxZ){
                        SimuNodes.at(i).Up = &(SimuNodes.at(i-BoxX*BoxY*(BoxZ-1)));
                    }
                }
//              /**/
            }
        }

/*************************************************************/

        std::vector<double> getMeshProp(int which, int index){
            std::vector<double> result(Num_Nodes,0.0);
#pragma parallel for
            for(int i = 0; i < Num_Nodes;i++){
                result.at(i) = (SimuNodes.at(i).getProp(which).at(index));
            }
            return result;
        }

        std::vector<double> getMeshCon(ELEMENT _element){
            std::vector<double> meshcon(Num_Nodes,0.0);
            for(int i = 0; i < Num_Nodes;i++){
                meshcon.at(i) = (SimuNodes.at(i).Con_Node.getCon(_element));
            }
            return meshcon;
        }

        int getDim(const int which){
            return Dimension.at(which);
        }

        double getNum_Nodes(){ // return the number of SimuNodes.at(i)s in mesh
            return SimuNodes.size();
        }

        double getStepLength(const int which){
            return StepLength.at(which);
        }

        int getNum_Prop(WHICHPARA which){
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


        std::vector<int> transCoord(int where){
            if(where<Num_Nodes){       
                std::vector<int>coord;
                coord.at(0) = where%BoxX;
                coord.at(1) = ((where-coord.at(0)) / BoxX)%BoxY;
                coord.at(2) = (where-(coord.at(0)+coord.at(1)*BoxX)/(BoxX*BoxY));
                return coord;
            }
            throw std::out_of_range("Not in mesh");
            return {};
        }

        int transCoord(std::vector<int> where){
            if(where.at(0)<BoxX&&where.at(1)<BoxY&&where.at(2)<BoxZ)
            return where.at(0)+where.at(1)*BoxX+where.at(2)*BoxX*BoxY;
            throw std::out_of_range("Not in mesh");
            return 0;
        }
        
        // find SimuNodes.at(i) in the mesh according to the coordinates
        MeshNode& findNode(int X, int Y, int Z){ 
            if(X<Dimension.at(0)&&Y<Dimension.at(1)&&Z<Dimension.at(2) && !(X<0) &&!(Y<0) &&!(Z<0)){
                return SimuNodes.at(X+Y*Dimension.at(0)+Z*Dimension.at(0)*Dimension.at(1));
            }
            else throw std::invalid_argument("Index Not in Mesh");
        }

        MeshNode& findNode(int AbsCoord){
            if(AbsCoord<Num_Nodes){
                return SimuNodes.at(AbsCoord);
            }
            else throw std::invalid_argument("Index Not in Mesh");
        }

        void updateCustomValue(std::vector<int> CustomValue){
            for(int i = 0; i < Num_Nodes; i++){
                SimuNodes.at(i).Custom_Value = CustomValue.at(i);
            }
        }

/*************************************************************/

        void updateMeshCon(std::vector<double> _val){
            for(int i = 0; i < Num_Nodes; i++){
                if(SimuNodes.at(i).getNum(WHICHPARA::CON) == 1)
                SimuNodes.at(i).Con_Node.updateEntry(_val.at(i));
                else throw std::invalid_argument("Exist more than one element");
            }
        }

        void updateMeshCon(ELEMENT _element,std::vector<double> _val){
            for(int i = 0; i < Num_Nodes; i++){
                SimuNodes.at(i).Con_Node.updateEntry(_element,_val.at(i));
            }
        }

        void updateNodeCon(int where, double _con){
            if(SimuNodes.at(where).getNum(WHICHPARA::CON) == 1)
            SimuNodes.at(where).Con_Node.updateEntry(_con);
            else throw std::invalid_argument("Exist more than one element");
        }
        
        void updateNodeCon(std::vector<int> where, double _con){
            updateNodeCon(transCoord(where) , _con);
        }

        void updateNodeCon(int where, ELEMENT _element, double _con){
            SimuNodes.at(where).Con_Node.updateEntry(_element,_con);
        } 

        void updateNodeCon(std::vector<int> where, ELEMENT _element, double _con){
            updateNodeCon(transCoord(where),_element,_con);
        }

/*************************************************************/
        void updateMeshPhs(int index, std::vector<double> _val){
            for(int i = 0; i < Num_Nodes; i++){
                SimuNodes.at(i).Phs_Node.updateEntry(index,_val.at(i));
            }
        }

        void updateNodePhs(int where,int index, double _phs){
            SimuNodes.at(where).Phs_Node.updateEntry(index,_phs);
        }

        void updateNodePhs(std::vector<int> where, int index, double _phs){
            updateNodePhs(transCoord(where),index,_phs);
        }



/*************************************************************/
        void showGlobalInfo(); // show the basic information of the mesh
        void showNodesProp(int which, int index); // show one of the properties of the SimuNodes.at(i)s inside the mesh
        void write_vtk_grid_values(int istep);
        void outFile(int istep);

/*************************************************************/

/**/    std::vector<double> Laplacian (int whichSTNCL, int whichPara){
            std::vector<double> result(Num_Nodes,0.0);
            double dx = getStepLength(WHICHDIM::X);
            double dy = getStepLength(WHICHDIM::Y);
            double dz = getStepLength(WHICHDIM::Z);
            if(whichSTNCL = STENCILE::FIVEPOINT){
            #pragma omp parallel for
                for(int i = 0; i < Num_Nodes; i++){
                    
                    double c = findNode(i).getProp(whichPara).at(0);
                    double f = findNode(i).Forward->getProp(whichPara).at(0);
                    double b = findNode(i).Backward->getProp(whichPara).at(0);
                    double l = findNode(i).Left->getProp(whichPara).at(0);
                    double r = findNode(i).Right->getProp(whichPara).at(0);
    
                    result.at(i) = ((f+b+l+r-4*c)/(dx*dy*dz));
                }
            }
            return result;
        }
};

void SimulationMesh::showGlobalInfo(){
    std::cout<<"SimulationMesh Properties:\n";
    std::cout<<"Mesh Size:\t\t"<<BoxX<<"\u0078"<<BoxY<<"\u0078"<<BoxZ<<"\n";
    std::cout<<"Number of Nodes:\t"<<Num_Nodes<<"\n";
    std::cout<<"Elements:\t\t";
    for(auto elements : SimuNodes.at(0).Con_Node.getElementList())
    std::cout<<elements<<" ";
    std::cout<<"\nNumber of Phase:\t"<<SimuNodes.at(0).getNum(WHICHPARA::PHSFRAC);

    std::cout<<"\n-----------------------------------------------------------------\n"<<std::endl;
}

void SimulationMesh::showNodesProp(int which, int index){ //which para, index of para
    if(which == WHICHPARA::CON)std::cout<<"Concentration of "<<SimuNodes.at(0).Con_Node.getElementList().at(index)<<"\n";
    if(which == WHICHPARA::PHSFRAC)std::cout<<"Order Parameter of "<<index<<" grain\n";
    if(which == WHICHPARA::CUSTOM)std::cout<<"Custom value \n";
    for(long i = 0; i < Num_Nodes; i++){
        std::cout<<std::fixed<<std::setprecision(10)<<SimuNodes.at(i).getProp(which).at(index)<<" ";
        if(i%BoxX==BoxX-1)std::cout<<"\n";
        if(i%(BoxX*BoxY) == BoxX*BoxY-1)std::cout<<"\n";
    }
    std::cout<<"-----------------------------------------------------------------\n"<<std::endl;
}

void SimulationMesh::outFile(int istep){

	char filename[128];
	sprintf(filename,"../output/Result_1/time_%04d.vtk", istep);

	std::ofstream outfile;
	outfile.open(filename);
	outfile<<"# vtk DataFile Version 2.0\n";
	outfile<<"time_10.vtk\n";
	outfile<<"ASCII\n";
	outfile<<"DATASET STRUCTURED_GRID\n";

	outfile<<"DIMENSIONS "<<BoxX<<"  "<<BoxY<<"  "<<BoxZ<<"\n";
	outfile<<"POINTS "<<Num_Nodes<<"   float\n";
	double dumx,dumy,dumz;
    #pragma parallel for collapse(3)
	for(int i=0;i<BoxX;i++)
		for(int j=0;j<BoxY;j++)
            for(int k = 0; k < BoxZ; k++)
		    {
		    	dumx = i*StepLength.at(0); dumy = j*StepLength.at(1); dumz = k*StepLength.at(2);
    
		    	outfile<<dumx<<"   "<<dumy<<"   "<<dumz<<"\n";
		    }

	outfile<<"POINT_DATA "<<Num_Nodes<<"\n";

    char varname[64];
    // int phs_num = getNum_Prop(WHICHPARA::PHSFRAC);
    // int con_num = getNum_Prop(WHICHPARA::CON);

    std::vector<std::vector<double>> meshphs;
    meshphs.reserve(getNum_Prop(WHICHPARA::PHSFRAC));
    for(int num = 0; num < getNum_Prop(WHICHPARA::PHSFRAC); num ++){
        meshphs.push_back(getMeshProp(WHICHPARA::PHSFRAC,num));
    }

    #pragma parallel for
    for(int num = 0; num < getNum_Prop(WHICHPARA::PHSFRAC); num++){
        sprintf(varname,"PHSFRAC_%01d",num);
	    outfile<<"SCALARS "<<varname<<"  float  1\n";
	    outfile<<"LOOKUP_TABLE default\n";
        #pragma parallel for
	    for(int i=0;i<Num_Nodes;i++)
	    	outfile<<meshphs.at(num).at(i)<<"\n";
    }

    std::vector<std::vector<double>> meshcon;
    for(int num = 0; num < getNum_Prop(WHICHPARA::CON); num ++){
        meshcon.push_back(getMeshProp(WHICHPARA::CON,num));
    }

    #pragma parallel for
	for(int num = 0; num < getNum_Prop(WHICHPARA::CON); num++){
        sprintf(varname,"CON_%01d",num);
	    outfile<<"SCALARS "<<varname<<"  float  1\n";
	    outfile<<"LOOKUP_TABLE default\n";

        #pragma parallel for
	    for(int i=0;i<Num_Nodes;i++)
	    	outfile<<meshcon.at(num).at(i)<<"\n";
    }

    meshphs.clear();
    meshcon.clear();

    outfile.close();

}

#endif