#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;


int main(){
    PFMTools::TimePoint start;
    start = PFMTools::now();

    std::string ProjName(toVTK_Path("../../TEST"));
    MeshNode node(std::vector<PhaseEntry>(69, Def_PhsEnt));
    SimulationMesh mesh({ 200, 200, 1 }, { 1, 1, 1 }, 0.01, node);
    // srand(time(0));

    // int R = 12,Xc,Yc;
    // for (int index = 0; index<mesh.getNum_Ent(WHICHPARA::PHSFRAC); index++){
    //     while (true){
    //         Xc = rand()%mesh.MeshX; Yc = rand()%mesh.MeshY;
    //         if (!(mesh.isOverlap(WHICHPARA::PHSFRAC, {Xc , Yc }, R, 3))){
    //             if (mesh.generateDisk(WHICHPARA::PHSFRAC, { Xc, Yc }, index, R)){
    //                 mesh.outCSV(ProjName,"cs5_64",Xc,Yc);
    //                 mesh.outVTKFilehead(ProjName, istep);
    //                 mesh.outVTKAve(ProjName, WHICHPARA::PHSFRAC, istep);
    //                 break;
    //             }
    //         }
    //     }
    // }
    int num = 69;
    vector<int> seeds = mesh.gnrtDiskSeeds(num,12,0);

    for(int i = 0; i < mesh.getNum_Ent(WHICHPARA::PHSFRAC); i++){
        mesh.generateDisk(WHICHPARA::PHSFRAC,{seeds.at(i),seeds.at(num+i)},i,12);
    }
    // mesh.generateDisk(WHICHPARA::PHSFRAC,{50,40},0,12);

    // vector<int>alphaX{ 12, 37, 62, 87, 112, 137, 162, 187, 12, 37, 62, 87, 112, 137, 162, 187, 12, 37, 62, 87, 112, 137, 162, 187, 12, 37, 62, 87, 112, 137, 162, 187, 12, 37, 62, 87, 112, 137, 162, 187, 12, 37, 62, 87, 112, 137, 162, 187, 12, 37, 62, 87, 112, 137, 162, 187, 12, 37, 62, 87, 112, 137, 162, 187 };
    // vector<int>alphaY{ 12, 12, 12, 12, 12, 12, 12, 12, 37, 37, 37, 37, 37, 37, 37, 37, 62, 62, 62, 62, 62, 62, 62, 62, 87, 87, 87, 87, 87, 87, 87, 87, 112, 112, 112, 112, 112, 112, 112, 112, 137, 137, 137, 137, 137, 137, 137, 137, 162, 162, 162, 162, 162, 162, 162, 162, 187, 187, 187, 187, 187, 187, 187, 187 };
    // vector<int>betaX{ 23, 46, 69, 92, 115, 138, 161, 184 };
    // vector<int>betaY{ 25, 50, 75, 100, 125, 150, 175 };

    // for (int index = 0; index<mesh.getNum_Ent(WHICHPARA::PHSFRAC); index++){
    //     mesh.generateDisk(WHICHPARA::PHSFRAC,{alphaX.at(index),alphaY.at(index)},index,12);
    // }

    int istep = 5;
    mesh.outVTKFilehead(ProjName, istep);
    mesh.outVTKAve(ProjName, WHICHPARA::PHSFRAC, istep);
    mesh.outVTKAll(ProjName, WHICHPARA::PHSFRAC, istep);

    PFMTools::RunTimeCounter(start);
    // auto duration = std::chrono::duration_cast< std::chrono::microseconds > (stop-start);
    // cout<<"\nTime taken by programme: "<<( double )duration.count()/1e6<<" seconds"<<endl;
    // cout<<"\nTime taken by programme: "<<( double )duration.count()/1e6<<" seconds"<<endl;

    // system("pause");
    return 0;
}

