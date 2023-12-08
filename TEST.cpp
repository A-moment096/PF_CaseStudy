#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;


int main(){
    PFMTools::TimePoint start = PFMTools::now(), dur = PFMTools::now();
    MeshNode _node(PhaseNode(std::vector<PhaseEntry>(1, Def_PhsEnt)));
    SimulationMesh mesh({ 200, 200, 1 }, { 1, 1, 1 }, 5.0e-3, _node);
    mesh.addEntry(WHICHPARA::CUSTOM, 1);
    double R = 12;
    mesh.addEntry(WHICHPARA::PHSFRAC, 68);
    vector<int> datas(PFMTools::readCSV("../TEST/DiskSeeds_69.csv"));
    for (int index = 0; index<mesh.getNum_Ent(WHICHPARA::PHSFRAC); index++){
        mesh.generateDisk(WHICHPARA::PHSFRAC, { datas.at(2*index), datas.at(2*index+1) }, index, R);
    }
    int nstep = 10000, nprint = 50;
    std::string _path(toVTK_Path("../../TEST/CS5_69Cell_TEST"));
cout<<"ready"<<endl;
    for (int istep = 0; istep<=nstep; istep++){
    #pragma omp parallel for
        for (auto &node:mesh.SimuNodes){
        #pragma omp critical
            {
                double sum2 = node.Phs_Node.sumPhsFrac2();
                node.Cust_Node.updateVal(0, node.Phs_Node.sumPhsFrac()*sum2-node.Phs_Node.sumPhsFrac3());
                node.Cust_Node.updateVal(1, sum2);
            }
            // double sum2 = node.Phs_Node.sumPhsFrac2();
            // node.Cust_Node.updateVal(0,node.Phs_Node.sumPhsFrac()*sum2-node.Phs_Node.sumPhsFrac3());
            // node.Cust_Node.updateVal(1,sum2);
        }


        if (istep%nprint==0){
            // mesh.outVTKFilehead(_path, istep);
            // mesh.outVTKWgtd(_path, WHICHPARA::PHSFRAC, istep);
            // mesh.outVTKAll(_path,WHICHPARA::PHSFRAC,istep);
            cout<<"Done Step: "<<istep;
            PFMTools::RunTimeCounter(dur, true);
            dur = PFMTools::now();

        }
    }

    PFMTools::RunTimeCounter(start, true);

    return 0;
}
