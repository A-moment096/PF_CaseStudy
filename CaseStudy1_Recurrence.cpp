#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

double dfdcon(double c){
    return 2*(c*(1-c)*(1-c)-c*c*(1-c));
}

int main()
{
    PhaseNode ConNode; // Node Properties
    PhaseSimulationMesh Box(ConNode);

    for(int i = 0; i<Box.getNumNodes(); i++){
        Box.findNode(i).updateNode(WHICHPARA::CON,0,0.4+0.01-double(rand()%200)/10000);
    }
    int nprint = 100;
    
    for(int istep = 0; istep<20000;istep++){
        Box.updateCustomValue(Box.Laplacian(STENCILE::FIVEPOINT,WHICHPARA::CON));

        for(int i = 0; i<Box.getNumNodes(); i++){
            double custom = Box.findNode(i).getProperties(WHICHPARA::CUSTOM).at(0);
            double cencon = Box.findNode(i).getProperties(WHICHPARA::CON).at(0);
            Box.findNode(i).updateNode(WHICHPARA::CUSTOM,0,dfdcon(cencon)-0.5*custom);
        }

        vector<double> temp;
        temp = Box.Laplacian(STENCILE::FIVEPOINT,WHICHPARA::CUSTOM);
        
        for(int i = 0; i<Box.getNumNodes(); i++){
            double cencon = Box.findNode(i).getProperties(WHICHPARA::CON).at(0);
            cencon+= 0.01*temp.at(i);
            if(cencon>=0.9999)cencon=0.9999;
            if(cencon<=0.0001)cencon=0.0001;
            Box.findNode(i).updateNode(WHICHPARA::CON,0,cencon);
        }
    
        vector<double> c(Box.getNumNodes(),0.4);
        for(int i = 0; i<Box.getNumNodes(); i++){
            c.at(i)= Box.findNode(i).getProperties(WHICHPARA::CON).at(0);
        }

        if(fmod(istep,nprint)==0)
        {   
            cout<<"Done Step: "<<istep<<endl;
            write_vtk_grid_values(Box.getDim(WHICHDIM::X),Box.getDim(WHICHDIM::Y),Box.getStepLength(WHICHDIM::X),Box.getStepLength(WHICHDIM::Y),istep,c);
        }  

    }
    return 0;
}