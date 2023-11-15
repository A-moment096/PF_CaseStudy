#pragma once
#ifndef PHASE_SIMULATION_MESH
#define PHASE_SIMULATION_MESH

#include <iostream>
#include <vector>
#include "MeshNode.hh"

enum WHICHDIM{X,Y,Z};
enum BOUNDCOND{PERIODIC,CONST,ADIABATIC};
enum STENCILE{FIVEPOINT=5,NINEPOINT=9};
enum WHICHPARA {CON,GRAIN,CUSTOM=99};


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

        std::vector<std::vector<double>> MeshCon;
        std::vector<std::vector<double>> MeshOrdPara;
        std::vector<double> MeshCstmVal;

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
            for(double i = 0; i < getNumNodes(); i++ ){
                SimuNodes.at(i).Forward = (i-BoxY < 0? &(SimuNodes.at(i)): &(SimuNodes.at(i-BoxY)) );
                SimuNodes.at(i).Backward = (i+BoxY >= getNumNodes()? &(SimuNodes.at(i)) : &(SimuNodes.at(i+BoxY)));
                SimuNodes.at(i).Left =(i-1 < 0? &(SimuNodes.at(i)) : &(SimuNodes.at(i-1)));
                SimuNodes.at(i).Right = (i+1 >= getNumNodes()? &(SimuNodes.at(i)) : &(SimuNodes.at(i+1)));
                SimuNodes.at(i).Down= (i-BoxX*BoxY < 0? &(SimuNodes.at(i)): &(SimuNodes.at(i-BoxX*BoxY)));
                SimuNodes.at(i).Up= (i+BoxX*BoxY >= getNumNodes()? &(SimuNodes.at(i)): &(SimuNodes.at(i+BoxX*BoxY)));

                if(whichBOUNDCOND == BOUNDCOND::PERIODIC){
                    if((int)i%(int)BoxX == 0){
                        SimuNodes.at(i).Left = &(SimuNodes.at(i+(BoxX-1)));
                    }
                    if((int)i%(int)BoxX == (BoxX-1)){
                        SimuNodes.at(i).Right = &(SimuNodes.at(i-(BoxX-1)));
                    }
                    if(((int)i%(int)(BoxX*BoxY))>=0 && ((int)i%(int)(BoxX*BoxY))<BoxX){
                        SimuNodes.at(i).Forward = &(SimuNodes.at(i+BoxX*(BoxY-1)));
                    }
                    if(((int)i%(int)(BoxX*BoxY))>=BoxX*(BoxY-1) && ((int)i%(int)(BoxX*BoxY))<BoxX*BoxY){
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
        // take properties of node to mesh, seperated by index
        void updatePropMesh(unsigned which){
            if(which == WHICHPARA::CON){
                for(unsigned i = 0; i < SimuNodes.at(0).Con_Node.getNums(); i++){
                    MeshCon.push_back({});
                    for(unsigned j = 0; j < getNumNodes(); j++){
                        MeshCon.at(i).push_back(SimuNodes.at(j).getProperties(which).at(i));
                    }
                }
            }
            if(which == WHICHPARA::GRAIN){
                for(unsigned i = 0; i < SimuNodes.at(0).getNums(which); i++){
                    MeshOrdPara.push_back({});
                    for(unsigned j = 0; j < getNumNodes(); j++){
                        MeshOrdPara.at(i).push_back(SimuNodes.at(j).getProperties(which).at(i));
                    }
                }
            }
            if(which == WHICHPARA::CUSTOM){
                for(unsigned i = 0; i < getNumNodes(); i++){
                    MeshCstmVal.push_back(SimuNodes.at(i).getProperties(WHICHPARA::CUSTOM).at(0));
                }
            }
        }

        void fillCstmMesh(){
            for(auto node : SimuNodes){
                MeshCstmVal.push_back(node.getProperties(WHICHPARA::CUSTOM).at(0));
            }
        }

        std::vector<double> getMeshProp(unsigned which, unsigned index){
            std::vector<double> meshcon(getNumNodes(),0.0);
            for(unsigned i = 0; i < getNumNodes();i++){
                meshcon.at(i) = (SimuNodes.at(i).getProperties(which).at(index));
            }
            return meshcon;
        }

        double getDim(const double which){
            return Dimension.at(which);
        }

        double getNumNodes(){ // return the number of SimuNodes.at(i)s in mesh
            return SimuNodes.size();
        }

        double getStepLength(const int which){
            return StepLength.at(which);
        }
        
        MeshNode& findNode(double X, double Y, double Z){ // find SimuNodes.at(i) in the mesh according to the coordinates
            if(X<Dimension.at(0)&&Y<Dimension.at(1)&&Z<Dimension.at(2) && !(X<0) &&!(Y<0) &&!(Z<0)){
                return SimuNodes.at(X+Y*Dimension.at(0)+Z*Dimension.at(0)*Dimension.at(1));
            }
            else throw std::invalid_argument("Index Not in Mesh");
        }

        MeshNode& findNode(double AbsCoord){
            if(AbsCoord<getNumNodes()){
                return SimuNodes.at(AbsCoord);
            }
            else throw std::invalid_argument("Index Not in Mesh");
        }

        void updateCustomValue(std::vector<double> CustomValue){
            for(int i = 0; i < getNumNodes(); i++){
                SimuNodes.at(i).updateNode(WHICHPARA::CUSTOM,0,CustomValue.at(i));
            }
        }

        void updateNodeProp(unsigned whichpara, unsigned index, std::vector<double> valueVec){
            for(int i = 0; i < getNumNodes(); i++){
                SimuNodes.at(i).updateNode(whichpara,index,valueVec.at(i));
            }
        }

        void showMeshInfo(); // show the basic information of the mesh
        void showNodesProp(unsigned which, unsigned index); // show one of the properties of the SimuNodes.at(i)s inside the mesh

/*************************************************************/

        std::vector<double> Laplacian (int whichSTNCL, int whichPara){
            std::vector<double> result(getNumNodes(),0.0);
            double dx = getStepLength(WHICHDIM::X);
            double dy = getStepLength(WHICHDIM::Y);
            double dz = getStepLength(WHICHDIM::Z);
            if(whichSTNCL = STENCILE::FIVEPOINT){
            #pragma omp parallel for
                for(int i = 0; i < getNumNodes(); i++){
                    
                    double c = findNode(i).getProperties(whichPara).at(0);
                    double f = findNode(i).Forward->getProperties(whichPara).at(0);
                    double b = findNode(i).Backward->getProperties(whichPara).at(0);
                    double l = findNode(i).Left->getProperties(whichPara).at(0);
                    double r = findNode(i).Right->getProperties(whichPara).at(0);
    
                    result.at(i) = ((f+b+l+r-4*c)/(dx*dy*dz));
                }
            }
            return result;
        }
};

void SimulationMesh::showMeshInfo(){
    std::cout<<"SimulationMesh Properties:"<<std::endl;
    std::cout<<"Mesh Size:\t\t"<<BoxX<<"\u0078"<<BoxY<<"\u0078"<<BoxZ<<std::endl;
    std::cout<<"Number of Nodes:\t"<<getNumNodes()<<std::endl;
}

void SimulationMesh::showNodesProp(unsigned which, unsigned index){ //which para, index of para
    if(which == WHICHPARA::CON)std::cout<<"Concentration of "<<index<<" element\n";
    if(which == WHICHPARA::GRAIN)std::cout<<"Order Parameter of "<<index<<" grain\n";
    if(which == WHICHPARA::CUSTOM)std::cout<<"Custom value \n";
    for(long i = 0; i < getNumNodes(); i++){
        std::cout<<std::fixed<<std::setprecision(10)<<SimuNodes.at(i).getProperties(which).at(index)<<" ";
        if(i%(int)BoxX==BoxX-1)std::cout<<"\n";
        if(i%(int)(BoxX*BoxY) == BoxX*BoxY-1)std::cout<<"\n";
    }
    std::cout<<std::endl;
}

#endif