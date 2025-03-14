#include "../include/PFM.hh"

using std::cout;
using std::endl;
using std::vector;

const double &&modulus = 5.0, &&lambda = 7.0, &&kappa = 60.0, &&mu = 40.0, &&xi = 1.5e3;

const double &&pi = 4.0*atan(1), &&GlobalTmp0 = 30/(lambda*lambda), &&GlobalTmp1 = 60*kappa/(lambda*lambda*xi);

double vn(double IntgPhi, double v_nA){
    return v_nA+GlobalTmp1*IntgPhi;
}

double term2(double phi, double sum){
    return -GlobalTmp0*(modulus*phi*(1-phi)*(1-2*phi)+2*kappa*sum);
}

int main(){
    PFMTools::TimePoint start, dur;

/*******************************************************************************************************/
    //Preparation
    //about file path, constants, parameters, mesh and nodes
    std::string _path(toVTK_Path("../../CS5_MultiCell_TEST"));

    MeshNode node(PhaseNode(std::vector<PhaseEntry>(1, Def_PhsEnt)));
    SimulationMesh mesh({ 200, 200, 1 }, { 1, 1, 1 }, 5.0e-3 ,node);

    const double &&v_nA = 0.5, &&R = 25;
    const double &&pi = 4.0*atan(1), &&Tmp0 = -2.0*mu/(pi*R*R);
    const int &&nstep = 3000, &&nprint = 50;
    const std::vector<double> v_0_1_A{0.5,-0.5};

    int simuflag = 1;
    if(simuflag ==1){
        mesh.addEntry(WHICHPARA::PHSFRAC,1);
    }
    else if(simuflag ==2){
        
    }
    

    if (simuflag==0||simuflag==1){
        double &&xc1 = mesh.MeshX/2-1.25*R;
        double &&yc1 = mesh.MeshY/2;
        double &&xc2 = mesh.MeshX/2+1.25*R;
        double &&yc2 = mesh.MeshY/2;

        for (int i = 0; i<mesh.MeshX; i++){
            for (int j = 0; j<mesh.MeshY; j++){
                if ((i-xc1)*(i-xc1)+(j-yc1)*(j-yc1)<=R*R){
                    mesh.updateNodeVal(WHICHPARA::PHSFRAC, { i, j, 0 }, 0, 0.999999);
                }
                if (simuflag==1&&(i-xc2)*(i-xc2)+(j-yc2)*(j-yc2)<=R*R){
                    mesh.updateNodeVal(WHICHPARA::PHSFRAC, { i, j, 0 }, 1, 0.999999);
                }
            }
        }
    }

    for (int istep = 0; istep<=nstep; istep++){
        if(istep == 500){mesh.setTimeStep(1.0e-2);}
        mesh.Laplacian(WHICHPARA::PHSFRAC);
        mesh.Gradient(WHICHPARA::PHSFRAC);

        for (int index = 0; index<mesh.getNum_Ent(WHICHPARA::PHSFRAC); index++){
            
            double Intg = 0.0, Intx = 0.0, Inty = 0.0;
            double VnX = 0.0, VnY = 0.0;

            for (auto &node0:mesh.SimuNodes){
                // Integral in cust
                double phi = node0.getVal(WHICHPARA::PHSFRAC, index);
                Intg += phi*phi;
                if(simuflag==2){
                    phi *= (node0.sumPhsFrac2()-phi*phi);
                    Intx += phi*node0.getGrad(WHICHPARA::PHSFRAC, index, DIM::DimX);
                    Inty += phi*node0.getGrad(WHICHPARA::PHSFRAC, index, DIM::DimY);
                }
            }

            double INTG = std::move(Intg);

            if(simuflag==2){
                VnX = std::move(vn(Intx, v_nA));
                VnY = std::move(vn(Inty, v_nA));
            }
            else if(simuflag==0||simuflag==1){
                VnX = v_0_1_A.at(index);
                VnY = 0;
            }

            for (auto &node1:mesh.SimuNodes){
                double &&Term1 = modulus*node1.getLap(WHICHPARA::PHSFRAC, index);
                double &&Term2 = term2(node1.getVal(WHICHPARA::PHSFRAC, index), node1.sumPhsFrac()*node1.sumPhsFrac2()-node1.sumPhsFrac3());
                double &&Term3 = Tmp0*node1.getVal(WHICHPARA::PHSFRAC, index)*(INTG-pi*R*R);
                double &&Term4 = (VnX*node1.getGrad(WHICHPARA::PHSFRAC, index, DIM::DimX)+VnY*node1.getGrad(WHICHPARA::PHSFRAC, index, DIM::DimY));
                node1.Phs_Node.updateDVal(index, Term1+Term2+Term3-Term4);
            }
        }

        mesh.iterateVal(WHICHPARA::PHSFRAC);

        if (istep%nprint==0){
            mesh.outVTKFilehead(_path, istep);
            mesh.outVTKAve(_path, WHICHPARA::PHSFRAC, istep);
            cout<<"Done Step: "<<istep;
            PFMTools::RunTimeCounter(dur);
            dur = PFMTools::now();
        }
    }

    PFMTools::RunTimeCounter(start);

//?range(boundary) check

//?sum_all function



    return 0;

}