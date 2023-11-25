#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

double dfdcon(double c){
    return 2*(c*(1-c)*(1-c)-c*c*(1-c));
}

int main(){
    MeshNode ConNode; // Node Properties
    SimulationMesh Box(ConNode);

    for (auto &node : Box.SimuNodes){
        node.Con_Node.updateVal(0, 0.4+0.01-double(rand()%200)/10000);
    }
    int nprint = 100;

    for (int istep = 0; istep<20001;istep++){
        Box.Laplacian(STENCILE::FIVEPOINT, WHICHPARA::CON);

    #pragma omp parallel for
        for (auto &node : Box.SimuNodes){
            double custom = node.Con_Node.getLap().at(0);
            double cencon = node.getProp(WHICHPARA::CON).at(0);
            node.Cust_Node.CustVal.at(0) = (dfdcon(cencon)-0.5*custom);
        }

        Box.Laplacian(STENCILE::FIVEPOINT, WHICHPARA::CUSTOM);

    #pragma omp parallel for
        for (auto &node : Box.SimuNodes){
            double cencon = node.getProp(WHICHPARA::CON).at(0);
            cencon += 0.01*node.Cust_Node.CustLap.at(0);
            Box.threshold(cencon, 0.0001, 0.9999);
            node.Con_Node.updateVal(0, cencon);
        }

        if (fmod(istep, nprint) == 0){
            cout<<"Done Step: "<<istep<<endl;
            Box.outFile(istep);
        }
    }
    return 0;
}