#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

int main(){
    auto start =  std::chrono::high_resolution_clock:: now();

    MeshNode node(PhaseNode (std::vector<PhaseEntry> (2,Def_PhsEnt)),Def_ConNode);
    SimulationMesh mesh({10,10,1},{0.5,0.5,1},node);

    srand(11111);
    for(int i = 0; i<mesh.Num_Nodes; i++){
        mesh.findNode(i).Con_Node.updateEntry(0,0.4+0.01-double(rand()%200)/10000);
    }

    mesh.Laplacian(STENCILE::FIVEPOINT,WHICHPARA::CON);
    mesh.showNodesProp(WHICHPARA::CON,0);
    
    auto stop =  std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start);
    cout<<"\nTime taken by programme: "<<(double) duration.count() / 1e6 << " seconds"<<endl;

    // system("pause");
    return 0;
}