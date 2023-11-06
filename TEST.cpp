#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

int main(){
    PhaseSimulationMesh Box(std::vector<long> {16,3,2});
    vector<PhaseNode> nodes(3);
    SpinodalSolver Solver;
    PhaseDomain NewDomain(nodes,Box);
    NewDomain.showBaseMeshProp();
    NewDomain.showGrainProp();
    Solver.showBaseMeshProp();
    Solver.showGrainProp();
    return 0;
}