#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

int main(){
    auto start =  std::chrono::high_resolution_clock:: now();
        
    PhaseNode phs_node(std::vector<PhaseEntry> (2,Def_PhsEnt));
    ConNode con_node;
    MeshNode node(phs_node,con_node);
    SimulationMesh mesh({100,100,1},{0.5,0.5,0},node);

    phs_node.addEntry(10);
    cout<<phs_node.getNums();

    // mesh.showGlobalInfo();
    // mesh.addEntry(WHICHPARA::PHSFRAC,10);
    // mesh.showGlobalInfo();

    auto stop =  std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start);
    cout<<"Time taken by programme: "<<(double) duration.count() / 1e6 << " seconds"<<endl;

    // system("pause");
    return 0;
}