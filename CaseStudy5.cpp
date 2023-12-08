#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

const double &&modulus = 5.0, &&lambda = 7.0, &&kappa = 60.0, &&mu = 40.0, &&xi = 1.5e3;

const double &&pi = 4.0*atan(1), &&GlobalTmp0 = 30/(lambda*lambda), &&GlobalTmp1 = 60*kappa/(lambda*lambda*xi);

double vn(double IntgPhi, double v_nA){
    return v_nA+GlobalTmp1*IntgPhi;
}

double term2_a( double phi, double sum){
    return -GlobalTmp0*(modulus*phi*(1-phi)*(1-2*phi)+2*kappa*sum);
}

double term2_b(double phi, double sum){
    return -GlobalTmp0*(0.5*modulus*phi*(1-phi)*(1-2*phi)+2*kappa*sum);
}

int main(){
    PFMTools::TimePoint start = PFMTools::now(), dur = PFMTools::now();

/*******************************************************************************************************/
    //Preparation
    //about file path, constants, parameters, mesh and nodes
    std::string _path(toVTK_Path("../../CS5_69Cell_TEST"));

    MeshNode _node(PhaseNode(std::vector<PhaseEntry>(1, Def_PhsEnt)));
    SimulationMesh mesh({ 200, 200, 1 }, { 1, 1, 1 }, 5.0e-3 ,_node);
    mesh.addEntry(WHICHPARA::CUSTOM,1);

    double R = 0.0;
    std::vector<double> Velo;

    int simuflag = 2;
    if(simuflag ==1){
        R = 25;
        mesh.addEntry(WHICHPARA::PHSFRAC,1);
        Velo.push_back(0.5);
        Velo.push_back(-0.5);
    }
    else if(simuflag ==2){
        R = 12;
        mesh.addEntry(WHICHPARA::PHSFRAC,68);
        vector<int> datas (PFMTools::readCSV("../TEST/DiskSeeds_69.csv"));
        for(int index = 0; index < mesh.getNum_Ent(WHICHPARA::PHSFRAC); index++){
            mesh.generateDisk(WHICHPARA::PHSFRAC,{datas.at(2*index),datas.at(2*index+1)},index,R);
            Velo.reserve(mesh.getNum_Ent(WHICHPARA::PHSFRAC));

            if(rand()%2==0){
                Velo.push_back(0.2);
            }
            else{
                Velo.push_back(-0.2);
            }
            if(index<5){
                for(auto &node : mesh.SimuNodes){
                    node.Phs_Node.Entrys.at(index).Weight = 1.3;
                }
            }
        }      
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

    const double &&pi = 4.0*atan(1), &&Tmp0 = -2.0*mu/(pi*R*R);
    const int &&nstep = 3000, &&nprint = 50;

    for (int istep = 0; istep<=nstep; istep++){
        if(istep == 500){mesh.setTimeStep(1.0e-2);}
        mesh.Laplacian(WHICHPARA::PHSFRAC);
        mesh.Gradient(WHICHPARA::PHSFRAC);
        #pragma omp parallel for
        for(auto &node : mesh.SimuNodes){
            #pragma omp critical
            {
                double sum2 = node.Phs_Node.sumPhsFrac2();
                node.Cust_Node.updateVal(0,node.Phs_Node.sumPhsFrac()*sum2-node.Phs_Node.sumPhsFrac3());
                node.Cust_Node.updateVal(1,sum2);
            }
        }
        for (int index = 0; index<mesh.getNum_Ent(WHICHPARA::PHSFRAC); index++){
            
            double Intg = 0.0, Intx = 0.0, Inty = 0.0;
            double VnX = 0.0, VnY = 0.0;

            for (auto &node0:mesh.SimuNodes){
                // Integral in cust
                double phi = node0.getVal(WHICHPARA::PHSFRAC, index);
                Intg += phi*phi;
                if(istep>50){
                    phi *= (node0.Cust_Node.getVal(1)-phi*phi);
                    Intx += phi*node0.getGrad(WHICHPARA::PHSFRAC, index, DIM::DimX);
                    Inty += phi*node0.getGrad(WHICHPARA::PHSFRAC, index, DIM::DimY);
                }
            }

            double INTG = std::move(Intg);
            if(istep >50){
                {
                    VnX = vn(Intx,Velo.at(index));
                    VnY = vn(Inty,Velo.at(index));
                    // VnY =0;
                }
            }

            for (auto &node1:mesh.SimuNodes){
                double &&Term1 = modulus*node1.getLap(WHICHPARA::PHSFRAC, index);
                double &&Term2 = index>=5?term2_a(node1.getVal(WHICHPARA::PHSFRAC, index),node1.getVal(WHICHPARA::CUSTOM,0))
                                         :term2_b(node1.getVal(WHICHPARA::PHSFRAC, index),node1.getVal(WHICHPARA::CUSTOM,0));
                double &&Term3 = Tmp0*node1.getVal(WHICHPARA::PHSFRAC, index)*(INTG-pi*R*R);
                double &&Term4 = (VnX*node1.getGrad(WHICHPARA::PHSFRAC, index, DIM::DimX)+VnY*node1.getGrad(WHICHPARA::PHSFRAC, index, DIM::DimY));
                node1.Phs_Node.updateDVal(index, Term1+Term2+Term3-Term4);
            }
        }

        mesh.iterateVal(WHICHPARA::PHSFRAC);

        if (istep%nprint==0){
            mesh.outVTKFilehead(_path, istep);
            mesh.outVTKWgtd(_path, WHICHPARA::PHSFRAC, istep);
            mesh.outVTKAll(_path,WHICHPARA::PHSFRAC,istep);
            cout<<"Done Step: "<<istep;
            PFMTools::RunTimeCounter(dur,true);
            dur = PFMTools::now();
        }
    }

    PFMTools::RunTimeCounter(start,true);

//?range(boundary) check

//?sum_all function



    return 0;

}