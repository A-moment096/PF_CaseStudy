#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

double dfdcon(double c){
    return 2*(c*(1-c)*(1-c)-c*c*(1-c));
}

int main()
{
    MeshNode ConNode; // Node Properties
    SimulationMesh Box(ConNode);

    for(int i = 0; i<Box.getNum_Nodes(); i++){
        Box.findNode(i).Con_Node.updateEntry(0.4+0.01-double(rand()%200)/10000);
    }
    int nprint = 100;
    
    for(int istep = 0; istep<20000;istep++){
        Box.updateCustomValue(Box.Laplacian(STENCILE::FIVEPOINT,WHICHPARA::CON));

        for(int i = 0; i<Box.getNum_Nodes(); i++){
            double custom = Box.findNode(i).getProp(WHICHPARA::CUSTOM).at(0);
            double cencon = Box.findNode(i).getProp(WHICHPARA::CON).at(0);
            Box.findNode(i).Custom_Value = (dfdcon(cencon)-0.5*custom);
        }

        vector<double> temp;
        temp = Box.Laplacian(STENCILE::FIVEPOINT,WHICHPARA::CUSTOM);
        
        for(int i = 0; i<Box.getNum_Nodes(); i++){
            double cencon = Box.findNode(i).getProp(WHICHPARA::CON).at(0);
            cencon+= 0.01*temp.at(i);
            if(cencon>=0.9999)cencon=0.9999;
            if(cencon<=0.0001)cencon=0.0001;
            Box.findNode(i).Con_Node.updateEntry(cencon);
        }
    
        vector<double> c(Box.getNum_Nodes(),0.4);
        for(int i = 0; i<Box.getNum_Nodes(); i++){
            c.at(i)= Box.findNode(i).getProp(WHICHPARA::CON).at(0);
        }

        if(fmod(istep,nprint)==0)
        {   
            cout<<"Done Step: "<<istep<<endl;
            write_vtk_grid_values0(Box.getDim(WHICHDIM::X),Box.getDim(WHICHDIM::Y),Box.getStepLength(WHICHDIM::X),Box.getStepLength(WHICHDIM::Y),istep,c);
        }  

    }
    return 0;
}