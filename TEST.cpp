#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

int main(){
    PhaseSimulationMesh TestMesh({8,8,1},DefaultNode);
    TestMesh.showMeshInfo();
    vector<double> customval(TestMesh.getNumNodes());
    for(int i = 0; i< TestMesh.getNumNodes(); i++){
        customval.at(i) = i;
    }
    // for(auto val : customval)cout<<val<<endl;

    TestMesh.updateNodeProp(WHICHPARA::CON,0,customval);
    TestMesh.updatePropMesh(WHICHPARA::CON);

    ConEntry NewEntry;
    cout<<NewEntry.getCon();
    cout<<NewEntry.getindex();
}