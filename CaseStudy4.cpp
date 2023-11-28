#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;


double dfdphi(double phi, double m){
    return phi*(1-phi)*(phi-0.5+m);
}
double dfdcon(double c, double sum3, double sum2){
    const double &&A = 16.0, &&B = 1.0;
    return B*(2*c+4*sum3-6*sum2)-2*A*c*(3*c-2*c*c-1);
}

double dfdeta(double c, double x, double sum2){
    const double &&B = 1.0;
    return 12*B*x*(-2*x+c*x+1-c+sum2);
}

int main(){
    auto start = std::chrono::high_resolution_clock::now();

    /*******************************************************************************************************/
        //Preparation
        //about file path, constants, parameters, mesh and nodes, the  
    std::string _path(toVTK_Path("../../NineGrainSint1"));

    MeshNode node(PhaseNode(std::vector<PhaseEntry>(2, Def_PhsEnt)), Def_ConNode);
    SimulationMesh mesh({ 300, 300, 1 }, { 0.03, 0.03, 1 }, node);

    double &&tau = 0.0003, &&epsilon_b = 0.01, &&mu = 1.0, &&kappa = 1.8;
    double &&delta = 0.02, &&aniso1 = 4.0, &&aniso2 = 6.0;
    double &&alpha = 0.9, &&gamma = 10.0, &&Teq = 1.0, &&theta = 0.2;
    
    int nstep = 4000, nprint = 50; double dtime = 1.0e-4;

    for(int istep = 0; istep < nstep; istep++){
        
    }




    RunTimeCounter(start);

    return 0;
}