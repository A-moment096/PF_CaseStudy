#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

int main(){
    PhaseNode phs_node(std::vector<PhaseEntry> {Def_PhsEnt,Def_PhsEnt});
    MeshNode node(phs_node);
    SimulationMesh mesh(node);
    vector<double> ones(mesh.getNum_Nodes(),1.0);
    vector<double> zeros(mesh.getNum_Nodes(),0.0);

    double cent = mesh.getDim(WHICHDIM::X)/2;
    double rad = 12;

    mesh.updateMeshPhs(0,ones);
    mesh.updateMeshPhs(1,zeros);

    for(double i = 0; i < mesh.getDim(WHICHDIM::X); i++){
        for(double j = 0; j < mesh.getDim(WHICHDIM::Y); j++){
            if( ((i-cent)*(i-cent)+(j-cent)*(j-cent)) < rad*rad ) {
                mesh.updateNodePhs({i,j,0},0,0.5);
                mesh.updateNodePhs({i,j,0},1,0.5);
            }
        }
    }
    std::vector<double> out = mesh.getUni_Prop(WHICHPARA::PHSFRAC);

    // mesh.outFile(0,mesh.getMeshProp(WHICHPARA::PHSFRAC,1));
    mesh.outFile(1,out);

    // for (int i = 0; i < mesh.getDim(WHICHDIM::X); i++)
    // {
    //     for (int j = 0; j < mesh.getDim(WHICHDIM::Y); j++)
    //     {
    //         if(i<1||i>=mesh.getDim(WHICHDIM::X)-1||j<1||j>=mesh.getDim(WHICHDIM::Y)-1){
    //             mesh.updateNodeCon({(double)i,(double)j,0},2.22222222222);
    //         }
    //     }
    // }
    // 
    // mesh.showNodesProp(WHICHPARA::PHSFRAC, 0);
    // mesh.showNodesProp(WHICHPARA::PHSFRAC, 1);
}