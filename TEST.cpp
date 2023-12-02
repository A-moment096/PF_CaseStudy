#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;


int main(){
    PFMTools::TimePoint start;
    start = PFMTools::now();


    std::string ProjName(toVTK_Path("../../TEST"));
    MeshNode node(std::vector<PhaseEntry>(16, Def_PhsEnt));
    SimulationMesh mesh({ 200, 200, 1 }, { 1, 1, 1 }, 0.01, node);

    vector<int> datas (PFMTools::readCSV("../TEST/DiskSeeds.csv"));
    // // int R = 12,Xc,Yc;
    // // for (int index = 0; index<mesh.getNum_Ent(WHICHPARA::PHSFRAC); index++){
    // //     while (true){
    // //         Xc = rand()%mesh.MeshX; Yc = rand()%mesh.MeshY;
    // //         if (!(mesh.isOverlap(WHICHPARA::PHSFRAC, {Xc , Yc }, R, 3))){
    // //             if (mesh.generateDisk(WHICHPARA::PHSFRAC, { Xc, Yc }, index, R)){
    // //                 mesh.outCSV(ProjName,"cs5_64",Xc,Yc);
    // //                 mesh.outVTKFilehead(ProjName, istep);
    // //                 mesh.outVTKAve(ProjName, WHICHPARA::PHSFRAC, istep);
    // //                 break;
    // //             }
    // //         }
    // //     }
    // // }
    // int num = 16;
    // vector<int> seeds = mesh.gnrtDiskSeeds(num,25,0);

    // std::size_t const half_size = seeds.size()/2;
    // std::vector<int> Xs(seeds.begin(),seeds.begin()+half_size);
    // std::vector<int> Ys(seeds.begin()+half_size,seeds.end());

    cout<<datas.size();

    for(int i = 0; i < mesh.getNum_Ent(WHICHPARA::PHSFRAC); i++){
        mesh.generateDisk(WHICHPARA::PHSFRAC,{datas.at(2*i),datas.at(2*i+1)},i,25);
    }

    int istep = 999;
    mesh.outVTKFilehead(ProjName, istep);
    mesh.outVTKAve(ProjName, WHICHPARA::PHSFRAC, istep);
    mesh.outVTKAll(ProjName, WHICHPARA::PHSFRAC, istep);
    // mesh.outCSV(ProjName,"DiskSeeds",Xs,Ys);

    PFMTools::RunTimeCounter(start,true);

    return 0;
}

