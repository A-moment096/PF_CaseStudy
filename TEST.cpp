#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;


int main(){
    PFMTools::TimePoint start;
    start = PFMTools::now();

    std::string ProjName(toVTK_Path("../../TEST"));int istep = 5;
    MeshNode node(std::vector<PhaseEntry>(69, Def_PhsEnt));
    SimulationMesh mesh({ 200,200, 1 }, { 1, 1, 1 }, 0.01, node);
    srand(time(0));

    int R = 12,Xc,Yc;
    for (int index = 0; index<mesh.getNum_Ent(WHICHPARA::PHSFRAC); index++){
        while (true){
            Xc = rand()%mesh.MeshX; Yc = rand()%mesh.MeshY;
            if (!(mesh.isOverlap(WHICHPARA::PHSFRAC, {Xc , Yc }, R, 3))){
                if (mesh.generateDisk(WHICHPARA::PHSFRAC, { Xc, Yc }, index, R)){
                    mesh.outCSV(ProjName,"cs5_64",Xc,Yc);
                    mesh.outVTKFilehead(ProjName, istep);
                    mesh.outVTKAve(ProjName, WHICHPARA::PHSFRAC, istep);
                    break;
                }
            }
        }
    }

    // mesh.outVTKFilehead(ProjName, istep);
    // mesh.outVTKAve(ProjName, WHICHPARA::PHSFRAC, istep);
    // mesh.outVTKAll(ProjName, WHICHPARA::PHSFRAC, istep);

    PFMTools::RunTimeCounter(start);
    // auto duration = std::chrono::duration_cast< std::chrono::microseconds > (stop-start);
    // cout<<"\nTime taken by programme: "<<( double )duration.count()/1e6<<" seconds"<<endl;
    // cout<<"\nTime taken by programme: "<<( double )duration.count()/1e6<<" seconds"<<endl;

    // system("pause");
    return 0;
}

