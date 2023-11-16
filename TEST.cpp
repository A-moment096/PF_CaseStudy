#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

int main(){
    ConNode con_node(std::vector<ELEMENT>{ELEMENT::Fe});
    MeshNode node(con_node);
    SimulationMesh mesh({8,8,1},node);
    for(unsigned i = 0; i < mesh.getNum_Nodes(); i++){
        mesh.updateNodeCon(i,0.5);
    }
    for (unsigned i = 0; i < mesh.getDim(WHICHDIM::X); i++)
    {
        for (unsigned j = 0; j < mesh.getDim(WHICHDIM::Y); j++)
        {
            if(i<1||i>=mesh.getDim(WHICHDIM::X)-1||j<1||j>=mesh.getDim(WHICHDIM::Y)-1){
                mesh.updateNodeCon({(double)i,(double)j,0},2.22222222222);
            }
        }
    }   
    mesh.showNodesProp(WHICHPARA::CON,0);
}