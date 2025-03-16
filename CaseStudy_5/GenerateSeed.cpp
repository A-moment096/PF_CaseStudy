#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;


int main(){
    PFMTools::TimePoint start;
    start = PFMTools::now();
    
    int num = 35;

    std::string ProjName(toVTK_Path("../../TEST"));
    MeshNode node(std::vector<PhaseEntry>(35, Def_PhsEnt));
    SimulationMesh mesh({ 200, 200, 1 }, { 1, 1, 1 }, 0.01, node);

    vector<int> seeds = mesh.gnrtDiskSeeds(num,18,0);

    std::size_t const half_size = seeds.size()/2;
    std::vector<int> Xs(seeds.begin(),seeds.begin()+half_size);
    std::vector<int> Ys(seeds.begin()+half_size,seeds.end());

    mesh.outCSV(ProjName,"35DskSds",Xs,Ys);

    vector<int> datas (PFMTools::readCSV("../TEST/35DskSds.csv"));

    cout<<datas.size();

    for(int i = 0; i < mesh.getNum_Ent(WHICHPARA::PHSFRAC); i++){
        mesh.generateDisk(WHICHPARA::PHSFRAC,{datas.at(2*i),datas.at(2*i+1)},i,18);
    }

    int istep = 999;
    mesh.outVTKFilehead(ProjName, istep);
    mesh.outVTKAve(ProjName, WHICHPARA::PHSFRAC, istep);
    mesh.outVTKAll(ProjName, WHICHPARA::PHSFRAC, istep);

    PFMTools::RunTimeCounter(start,true);

    return 0;
}

