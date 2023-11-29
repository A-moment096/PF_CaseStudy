#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

using Timer = std::chrono::high_resolution_clock ;



int main(){
    auto start = Timer::now();


    MeshNode node(Def_ConNode);
    SimulationMesh mesh({ 3,3,1}, { 0.5, 0.5, 1 },0.01, node);
    

    RunTimeCounter(start);
    // #pragma omp parallel for
        for (int i = 0; i<mesh.Num_Nodes; ++i){
            mesh.updateNodeVal(WHICHPARA::CON,i,0,i);
        }

    // mesh.showNodesProp(WHICHPARA::CON,0);

    // cout<<mesh(4).getNbhd(WHICHDIR::DirF)->Con_Node.getVal(0)<<endl;
    // cout<<mesh(4).getNbhd(WHICHDIR::DirB)->Con_Node.getVal(0)<<endl;
    // cout<<mesh(4).getNbhd(WHICHDIR::DirL)->Con_Node.getVal(0)<<endl;
    // cout<<mesh(4).getNbhd(WHICHDIR::DirR)->Con_Node.getVal(0)<<endl;
    // cout<<mesh(4).getNbhd(WHICHDIR::DirU)->Con_Node.getVal(0)<<endl;
    // cout<<mesh(4).getNbhd(WHICHDIR::DirD)->Con_Node.getVal(0)<<endl;

    // std::string ProjName("../../TEST");int istep = 5;
    // ProjName = toVTK_Path(ProjName);
    // // mesh.outVTK(ProjName, istep);
    // mesh.outCSV(ProjName,"testcsv",0,0);

    double  TESTDouble = 500;

    TESTDouble *= 0 - TESTDouble*TESTDouble;
    cout<<500*500*500;

    cout<<TESTDouble<<endl;
    RunTimeCounter(start);
    // auto duration = std::chrono::duration_cast< std::chrono::microseconds > (stop-start);
    // cout<<"\nTime taken by programme: "<<( double )duration.count()/1e6<<" seconds"<<endl;
    // cout<<"\nTime taken by programme: "<<( double )duration.count()/1e6<<" seconds"<<endl;

    // system("pause");
    return 0;
}

