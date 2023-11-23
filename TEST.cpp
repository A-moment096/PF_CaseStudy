#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

int main(){
    auto start =  std::chrono::high_resolution_clock:: now();

    MeshNode node(PhaseNode (std::vector<PhaseEntry> (2,Def_PhsEnt)),Def_ConNode);
    SimulationMesh mesh({10,10,1},{0.5,0.5,0},node);
    
    auto stop =  std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start);
    cout<<"Time taken by programme: "<<(double) duration.count() / 1e6 << " seconds"<<endl;

    // system("pause");
    return 0;
}