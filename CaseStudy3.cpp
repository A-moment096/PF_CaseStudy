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


    double Dvol = 0.040, Dvap = 0.002, Dsurf = 16.0, Dgb = 1.6, kappa_rho = 5.0, kappa_eta = 2.0, L = 10.0;
    int nstep = 5000, nprint = 50; double dtime = 0.0001;

// Initializer
int simuflag = 1;
    if(simuflag == 0){    
    double cent = mesh.getDim(WHICHDIM::X)/2;
    double rad1 = 20, rad2 = 10;
    int Px = mesh.getDim(0)/2, Py1 = 40, Py2 = 70;

        #pragma omp parallel for collapse(2)
            for(int i = 0; i < mesh.getDim(WHICHDIM::X); ++i){
                for(int j = 0; j < mesh.getDim(WHICHDIM::Y); j++){
                    double &&dis1 = ((i-Px)*(i-Px)+(j-Py1)*(j-Py1));
                    double &&dis2 = ((i-Px)*(i-Px)+(j-Py2)*(j-Py2));
                    if( dis1<=rad1*rad1 ) {
                        mesh.updateNodeCon({i,j,0},0.9999);
                        mesh.updateNodePhs({i,j,0},0,0.9999);
                    }
                    if(dis2<=rad2*rad2){
                        mesh.updateNodeCon({i,j,0},0.9999);
                        mesh.updateNodePhs({i,j,0},0,0);
                        mesh.updateNodePhs({i,j,0},1,0.9999);
                    }
                }
        }
    }
    if(simuflag == 1){
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
            for(int coord = 0; coord < mesh.getNum_Prop(WHICHPARA::PHSFRAC); ++coord)
                for(int i = 0; i < mesh.getDim(WHICHDIM::X); ++i)
                    for(int j = 0; j < mesh.getDim(WHICHDIM::Y); ++j){
                        double dis = (i-particleCoord.at(2*coord))*(i-particleCoord.at(2*coord))+(j-particleCoord.at(2*coord+1))*(j-particleCoord.at(2*coord+1));
                        if(coord<6){
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
    mesh.showGlobalInfo();
    }
    mesh.showGlobalInfo();

    mesh.outFile(1);    

    auto stop =  std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start);
    cout<<"Time taken by programme: "<<(double) duration.count() / 1e6 << " seconds"<<endl;

    // system("pause");
    return 0;
}