#include "../include/PFM.hh"

using std::cout;
using std::endl;
using std::vector;

const double &&tau = 0.0003, &&epsilon_b = 0.01, &&mu = 1.0, &&kappa = 1.8;
const double &&delta = 0.02, &&aniso = 4;
const double &&alpha = 0.9, &&GaMMA = 10.0, &&Teq = 1.0, &&theta0 = 0.2;
const double &&pi = 4*atan(1);

double dfdphi(double phi, double m){
    return phi*(1-phi)*(phi-0.5+m);
}

double epsilon(double theta){
    return epsilon_b*(1+delta*cos(aniso*(theta-theta0)));
}

double dedtheta(double theta){
    return epsilon_b*(aniso*delta*sin(aniso*(theta0-theta)));
}

double m_cal(double T){
    return alpha/pi*atan(GaMMA*(Teq-T));
}


int main(){
    auto start = std::chrono::high_resolution_clock::now();

    /*******************************************************************************************************/
        //Preparation
        //about file path, constants, parameters, mesh and nodes, the
    std::string _path(toVTK_Path("./"));

    MeshNode node;
    SimulationMesh mesh({ 300, 300, 1 }, { 0.030, 0.030, 1 },1.0e-4 ,std::move(node));
    mesh.addEntry(WHICHPARA::CUSTOM, 2);
    cout<<"ready"<<endl;

    int nstep = 4000, nprint = 50; double dtime = 1.0e-4;
    #pragma omp parallel for collapse(2)
    for (int i = 0; i<mesh.MeshX; i++){
        for (int j = 0; j<mesh.MeshY; j++){
            if ((i-mesh.MeshX/2)*(i-mesh.MeshX/2)+(j-mesh.MeshY/2)*(j-mesh.MeshY/2)<=5.0*5.0){
                mesh.updateNodeVal(WHICHPARA::PHSFRAC,{ i, j, 0 }, 0, 1.0);
            }
        }
    }

    cout<<"ready to go"<<endl;

    for (int istep = 0; istep<= nstep; istep++){
        mesh.Laplacian(WHICHPARA::TEMP);
        mesh.Laplacian(WHICHPARA::PHSFRAC);
        mesh.Gradient(WHICHPARA::PHSFRAC);
        
        #pragma omp parallel for
        for (auto &node:mesh.SimuNodes){
            double Theta = atan2(node.getGrad(WHICHPARA::PHSFRAC, 0, DIM::DimY),node.getGrad(WHICHPARA::PHSFRAC, 0, DIM::DimX));
            double Epsilon = epsilon(Theta);
            node.Cust_Node.updateVal(0, Epsilon);
            node.Cust_Node.updateVal(1, Epsilon*dedtheta(Theta)*node.Phs_Node.getGrad(0, DIM::DimX));
            node.Cust_Node.updateVal(2, Epsilon*dedtheta(Theta)*node.Phs_Node.getGrad(0, DIM::DimY));
        }
        mesh.GradientY(WHICHPARA::CUSTOM, 1);
        mesh.GradientX(WHICHPARA::CUSTOM, 2);

        double m = 0.0;
        #pragma omp parallel for
        for (auto &node:mesh.SimuNodes){
            m = m_cal(node.Temp_Node.getVal(0));
            double dummy1 = dfdphi(node.Phs_Node.getVal(0), std::move(m));
            double dummy2 = node.getGrad(WHICHPARA::CUSTOM, 1, DIM::DimY)-node.getGrad(WHICHPARA::CUSTOM, 2, DIM::DimX);
            double dummy3 = dummy1+dummy2+node.Cust_Node.getVal(0)*node.Cust_Node.getVal(0)*node.getLap(WHICHPARA::PHSFRAC, 0);
            dummy3 /= tau;
            double dummy4 = dummy3*dtime+node.Phs_Node.getVal(0);
            PFMTools::threshold(dummy4,0.00001,0.99999);
            node.Phs_Node.updateVal(0,dummy4);
            node.Temp_Node.updateVal(0, (dummy3*kappa+node.Temp_Node.getLap(0))*dtime+node.Temp_Node.getVal(0));
        }
        if (istep%nprint==0){
            mesh.outVTKFilehead(_path, istep);
            mesh.outVTKAll(_path, WHICHPARA::PHSFRAC, istep);
            mesh.outVTKAll(_path, WHICHPARA::TEMP, istep);
            cout<<"Done Step: "<<istep;
            PFMTools::RunTimeCounter(start);
        }
    }

    PFMTools::RunTimeCounter(start);

    return 0;
}