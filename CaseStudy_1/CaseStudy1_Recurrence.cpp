#include "../include/PFM.hh"

using std::cout;
using std::endl;
using std::vector;

double dfdcon(double c){
    return 2*(c*(1-c)*(1-c)-c*c*(1-c));
}

int main(){
    auto start( std::chrono::high_resolution_clock::now());
    MeshNode ConNode; // Node Properties
    SimulationMesh Box(ConNode);

    std::string path(toVTK_Path("../../CaseStudy1_Recurrence"));

    for (auto &node : Box.SimuNodes){
        node.Con_Node.updateVal(0, 0.4+0.01-double(rand()%200)/10000);
    }
    int nprint = 100;

    for (int istep = 0; istep<20001;istep++){
        Box.Laplacian(WHICHPARA::CON);

    #pragma omp parallel for
        for (auto &node : Box.SimuNodes){
            double custom = node.Con_Node.getLap(0);
            double cencon = node.getVal(WHICHPARA::CON,0);
            node.Cust_Node.updateVal(0,(dfdcon(cencon)-0.5*custom));
        }

        Box.Laplacian( WHICHPARA::CUSTOM);

    #pragma omp parallel for
        for (auto &node : Box.SimuNodes){
            double cencon = node.getVal(WHICHPARA::CON,0);
            cencon += 0.01*node.Cust_Node.getLap(0);
            Box.threshold(cencon, 0.0001, 0.9999);
            node.Con_Node.updateVal(0, cencon);
        }

        if (fmod(istep, nprint) == 0){
            cout<<"Done Step: "<<istep<<endl;
            Box.outVTK(path,istep);
        }
    }

    RunTimeCounter(start);
    return 0;
}