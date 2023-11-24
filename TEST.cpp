#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

int main(){
    auto start =  std::chrono::high_resolution_clock:: now();

    MeshNode node(PhaseNode (std::vector<PhaseEntry> (2,Def_PhsEnt)),Def_ConNode);
    SimulationMesh mesh({100,100,1},{0.5,0.5,1},node);

    {
        mesh.addEntry(WHICHPARA::PHSFRAC,7);
        vector<int> particleCoord (18,0);
        particleCoord.at(0) = 29;
        particleCoord.at(1) = 50;
        particleCoord.at(2) = 50;
        particleCoord.at(3) = 50;
        particleCoord.at(4) = 71;
        particleCoord.at(5) = 50;
        particleCoord.at(6) = 50;
        particleCoord.at(7) = 29;
        particleCoord.at(8) = 50;
        particleCoord.at(9) = 71;

        particleCoord.at(10) = 39;
        particleCoord.at(11) = 39;
        particleCoord.at(12) = 61;
        particleCoord.at(13) = 39;
        particleCoord.at(14) = 39;
        particleCoord.at(15) = 61;
        particleCoord.at(16) = 61;
        particleCoord.at(17) = 61;
        double rad1 = 10.0;
        double rad2 = 5.0;
        #pragma omp parallel for collapse(2)
            for(int coord = 0; coord < mesh.getNum_Ent(WHICHPARA::PHSFRAC); ++coord)
                for(int i = 0; i < mesh.MeshX; ++i)
                    for(int j = 0; j < mesh.MeshY; ++j){
                        double dis = (i-particleCoord.at(2*coord))*(i-particleCoord.at(2*coord))+(j-particleCoord.at(2*coord+1))*(j-particleCoord.at(2*coord+1));
                        if(coord<5){
                            if(dis<=rad1*rad1){
                                mesh.updateNodeCon({i,j,0},0.9999);
                                mesh.updateNodePhs({i,j,0},coord,0.9999);
                            }
                        }
                        else{
                            if(dis<=rad2*rad2){
                                mesh.updateNodeCon({i,j,0},0.9999);
                                mesh.updateNodePhs({i,j,0},coord,0.9999);
                            }
                        }
                    }
    }

    mesh.outFilehead(0);
    mesh.outAll(WHICHPARA::CON,0);
    mesh.outAll(WHICHPARA::PHSFRAC,0);
    mesh.outAve(WHICHPARA::PHSFRAC,0);

    auto stop =  std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start);
    cout<<"\nTime taken by programme: "<<(double) duration.count() / 1e6 << " seconds"<<endl;

    // system("pause");
    return 0;
}

