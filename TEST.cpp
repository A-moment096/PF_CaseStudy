#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

int main(){
    auto start =  std::chrono::high_resolution_clock:: now();
    
    PhaseNode phs_node
    (std::vector<PhaseEntry> {Def_PhsEnt,Def_PhsEnt});
    PhaseEntry ent;
    MeshNode node(phs_node);
    SimulationMesh mesh({16,16,16},node);
//     vector<double> ones(mesh.getNum_Nodes(),1.0);
//     vector<double> zeros(mesh.getNum_Nodes(),0.0);

// //     double cent = mesh.getDim(WHICHDIM::X)/2;
//     double rad = 1;

//     mesh.updateMeshPhs(0,ones);
//     mesh.updateMeshPhs(1,zeros);

// #pragma omp parallel for collapse(3)
//     for(double i = 0+cent-2*rad; i < mesh.getDim(WHICHDIM::X)-cent+2*rad; i++){
//         for(double j = 0+cent-2*rad; j < mesh.getDim(WHICHDIM::Y)-cent+2*rad; j++){
//             for(double k = 0+cent-2*rad; k < mesh.getDim(WHICHDIM::Z)-cent+2*rad; k++)
//             if( ((i-cent)*(i-cent)+(j-cent)*(j-cent)) < rad*rad ) {
//                 mesh.updateNodePhs({i,j,k},0,0.5);
//                 mesh.updateNodePhs({i,j,k},1,0.5);
//             }
//         }
//     }
    // std::vector<double> out = mesh.getUni_Prop(WHICHPARA::PHSFRAC);
    

    // mesh.outFile(0,mesh.getMeshProp(WHICHPARA::PHSFRAC,1));

    // for (int i = 0; i < mesh.getDim(WHICHDIM::X); i++)
    // {
    //     for (int j = 0; j < mesh.getDim(WHICHDIM::Y); j++)
    //     {
    //         if(i<1||i>=mesh.getDim(WHICHDIM::X)-1||j<1||j>=mesh.getDim(WHICHDIM::Y)-1){
    //             mesh.updateNodeCon({(double)i,(double)j,0},2.22222222222);
    //         }
    //     }
    // }
    // mesh.outFile(1);
    
    // mesh.showNodesProp(WHICHPARA::PHSFRAC, 0);
    // mesh.showNodesProp(WHICHPARA::PHSFRAC, 1);
    auto stop =  std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start);
    cout<<"Time taken by programme: "<<(double) duration.count() / 1e6 << " seconds"<<endl;

    return 0;
}