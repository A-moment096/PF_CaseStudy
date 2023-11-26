#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

using Timer = std::chrono::high_resolution_clock ;



int main(){
    auto start = Timer::now();


    MeshNode node(Def_ConNode);
    SimulationMesh mesh({ 100,100,1 }, { 0.5, 0.5, 1 }, node);
    
        // double rad1 = 20, rad2 = 10;
        // int Px = mesh.MeshX/2, Py1 = 40, Py2 = 70;

    RunTimeCounter(start);
    #pragma omp parallel for
        for (int i = 0; i<mesh.Num_Nodes; ++i){
            mesh.updateNodeCon(i,i);
        }
    vector<double> testdouble(10000,225);
    double ts = 100;
    
    cout<<mesh(0).getNbhd(WHICHDIR::DirF)->Con_Node.getVal(0)<<endl;
    cout<<mesh(0).getNbhd(WHICHDIR::DirB)->Con_Node.getVal(0)<<endl;
    cout<<mesh(0).getNbhd(WHICHDIR::DirL)->Con_Node.getVal(0)<<endl;
    cout<<mesh(0).getNbhd(WHICHDIR::DirR)->Con_Node.getVal(0)<<endl;
    cout<<mesh(0).getNbhd(WHICHDIR::DirU)->Con_Node.getVal(0)<<endl;
    cout<<mesh(0).getNbhd(WHICHDIR::DirD)->Con_Node.getVal(0)<<endl;
    cout<<sizeof(mesh)<<" "<<sizeof(testdouble)<<" "<<sizeof(ts);

    // std::string ProjName("../../TEST");int istep = 5;
    // ProjName = toVTK_Path(ProjName);
    // mesh.outVTK(ProjName, istep);

    RunTimeCounter(start);
    // auto duration = std::chrono::duration_cast< std::chrono::microseconds > (stop-start);
    // cout<<"\nTime taken by programme: "<<( double )duration.count()/1e6<<" seconds"<<endl;
    // cout<<"\nTime taken by programme: "<<( double )duration.count()/1e6<<" seconds"<<endl;

    // system("pause");
    return 0;
}

