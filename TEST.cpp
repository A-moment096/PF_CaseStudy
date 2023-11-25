#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

int main(){
    auto start = std::chrono::high_resolution_clock::now();


    MeshNode node(PhaseNode(std::vector<PhaseEntry>(2, Def_PhsEnt)), Def_ConNode);
    SimulationMesh mesh({ 100, 100, 1 }, { 0.5, 0.5, 1 }, node);
    {
        double rad1 = 20, rad2 = 10;
        int Px = mesh.MeshX/2, Py1 = 40, Py2 = 70;

    #pragma omp parallel for collapse(2)
        for (int i = 0; i<mesh.MeshX; ++i){
            for (int j = 0; j<mesh.MeshY; j++){
                double &&dis1 = ((i-Px)*(i-Px)+(j-Py1)*(j-Py1));
                double &&dis2 = ((i-Px)*(i-Px)+(j-Py2)*(j-Py2));
                if (dis1<=rad1*rad1){
                    mesh.updateNodeCon({ i, j, 0 }, 0.9999);
                    mesh.updateNodePhs({ i, j, 0 }, 0, 0.9999);
                }
                if (dis2<=rad2*rad2){
                    mesh.updateNodeCon({ i, j, 0 }, 0.9999);
                    mesh.updateNodePhs({ i, j, 0 }, 0, 0);
                    mesh.updateNodePhs({ i, j, 0 }, 1, 0);
                }
            }
        }
    }


    std::string ProjName("../../TEST");int istep = 5;
    ProjName = toVTK_Path(ProjName);
mesh.outFile(ProjName,istep);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast< std::chrono::microseconds > (stop-start);
    cout<<"\nTime taken by programme: "<<( double )duration.count()/1e6<<" seconds"<<endl;

    // system("pause");
    return 0;
}

