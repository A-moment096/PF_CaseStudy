#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

double dfdcon(double c){
    return 2*(c*(1-c)*(1-c)-c*c*(1-c));
}

int main()
{
    PhaseNode EmptyNode;
    PhaseSimulationMesh Box(EmptyNode);

    for(int i = 0; i<Box.getNumNodes(); i++){
        Box.findNode(i).updateNode(WHICHPARA::CON,0,0.4+0.01-double(rand()%200)/10000);
    }
    int nprint = 100;

    PhaseSolver Solver;
    
    for(int istep = 0; istep<20000;istep++){
        Box.updateCustomValue(Solver.Laplacian(STENCILE::FIVEPOINT,Box,WHICHPARA::CON));

        for(int i = 0; i<Box.getNumNodes(); i++){
            double custom = Box.findNode(i).getProperties(WHICHPARA::CUSTOM).at(0);
            double cencon = Box.findNode(i).getProperties(WHICHPARA::CON).at(0);
            Box.findNode(i).updateNode(WHICHPARA::CUSTOM,0,dfdcon(cencon)-0.5*custom);
        }
        // Box.updateCustomValue(Solver.Laplacian(STENCILE::FIVEPOINT,Box,WHICHPARA::CUSTOM));

        vector<double> temp;
        temp = Solver.Laplacian(STENCILE::FIVEPOINT,Box,WHICHPARA::CUSTOM);
        
        for(int i = 0; i<Box.getNumNodes(); i++){
            double cencon = Box.findNode(i).getProperties(WHICHPARA::CON).at(0);
            cencon+= 0.01*temp.at(i);
            if(cencon>=0.9999)cencon=0.9999;
            if(cencon<=0.0001)cencon=0.0001;
            Box.findNode(i).updateNode(WHICHPARA::CON,0,cencon);
        }
    
        vector<vector<double>> c(64, vector<double>(64,0.4));
        for(int i = 0; i<Box.getDim(WHICHDIM::X); i++){
            for(int j = 0; j<Box.getDim(WHICHDIM::Y); j++)
            c.at(i).at(j) = Box.findNode(i,j,0).getProperties(WHICHPARA::CON).at(0);
        }

        if(fmod(istep,nprint)==0)
        {   
            cout<<"Done Step: "<<istep<<endl;
            write_vtk_grid_values(Box.getDim(WHICHDIM::X),Box.getDim(WHICHDIM::Y),Box.getStepLength(WHICHDIM::X),Box.getStepLength(WHICHDIM::Y),istep,c);
        }  

    }
    return 0;
}