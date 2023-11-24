#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

double phi(double c){
    return c*c*c*(10-15*c+6*c*c);
}

int main(){
    auto start =  std::chrono::high_resolution_clock:: now();

    MeshNode node(PhaseNode (std::vector<PhaseEntry> (2,Def_PhsEnt)),Def_ConNode);
    SimulationMesh mesh({100,100,1},{0.5,0.5,1},node);

    double Dvol = 0.040, Dvap = 0.002, Dsurf = 16.0, Dgb = 1.6, kappa_rho = 5.0, kappa_eta = 2.0, L = 10.0;
    int nstep = 5000, nprint = 50; double dtime = 1e-4;
    double coefm = 5.0, coefk = 2.0, coefl = 5.0;

// Initializer
int simuflag = 0;
    if(simuflag == 0){    
    double cent = mesh.MeshX/2;
    double rad1 = 20, rad2 = 10;
    int Px = mesh.MeshX/2, Py1 = 40, Py2 = 70;

        #pragma omp parallel for collapse(2)
            for(int i = 0; i < mesh.MeshX; ++i){
                for(int j = 0; j < mesh.MeshY; j++){
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
    double cent = mesh.MeshX/2;
    double rad1 = 10;
    int Px = mesh.MeshX/2, Py = mesh.MeshY/2;

        #pragma omp parallel for collapse(2)
            for(int i = 0; i < mesh.MeshX; ++i){
                for(int j = 0; j < mesh.MeshY; j++){
                    double &&dis = ((i-Px)*(i-Px)+(j-Py)*(j-Py));
                    if( dis >= rad1*rad1 && i >= Px) {
                        mesh.updateNodeCon({i,j,0},0.9999);
                        mesh.updateNodePhs({i,j,0},0,0.9999);
                    }
                    if( dis >= rad1*rad1 && i <= Px) {
                        mesh.updateNodeCon({i,j,0},0.9999);
                        mesh.updateNodePhs({i,j,0},1,0.9999);
                    }
                }
        }
    }

    if(simuflag == 2){
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

    mesh.outFile(0);    

for (int istep = 0; istep < nstep; ++istep)
{
    mesh.Laplacian(WHICHPARA::CON);
    mesh.Laplacian(WHICHPARA::PHSFRAC);

    double A = 16.0, B = 1.0;
    for(auto &node : mesh.SimuNodes){
        double dfdcon, c = node.Con_Node.getCon().at(0); // // // // // node.getCon().at(0)
        dfdcon = B*(2*c+4*node.sumPhsFrac3() - 6*node.sumPhsFrac2())- 2*A*(3*c*c*-2*c*c*c-c);
        node.Custom_Value.at(0) = (dfdcon - 0.5*coefm* (node.getLap(WHICHPARA::CON).at(0)) );
    }

    for(auto &node: mesh.SimuNodes){
        double dfdeta, c = node.Con_Node.getCon().at(0);
        for(int i = 0; i < node.getNum_Ent(WHICHPARA::PHSFRAC); ++i){
            double x = node.Phs_Node.getPhsFrac().at(i);
            dfdeta = 12*B*(x*(-2*x+c*x+1-c+node.sumPhsFrac2()));
            double dummy = x-dtime*coefl*(dfdeta-0.5*coefk*node.Phs_Node.getLap().at(i));
            mesh.threshold(dummy,0.0001,0.9999);
            node.Phs_Node.updateEntry(i,dummy);
        }
    }

    mesh.Laplacian(WHICHPARA::CUSTOM);

    double Diffu = 0;
    
    for(auto &node : mesh.SimuNodes){
        double sum = node.sumPhsFrac()*node.sumPhsFrac() - node.sumPhsFrac2();
        double c = node.Con_Node.getCon().at(0);
        Diffu = Dvol*(phi(c))+Dvap*(1-phi(c))+Dsurf*c*(1-c)+Dgb*sum;

        double dumy = c + dtime*Diffu*node.CustLap;
        
        mesh.threshold(dumy,0.0001,0.9999);

        node.Con_Node.updateEntry(0,dumy);
    }
    
    if(fmod(istep,nprint)==0)
        {   
            mesh.outFile(istep);    
            cout<<"Done Step: "<<istep<<endl;
        }  
}



    auto stop =  std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start);
    cout<<"\nTime taken by programme: "<<(double) duration.count() / 1e6 << " seconds"<<endl;

    return 0;
}

