#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

<<<<<<< HEAD
int main()
{
    PhaseNode EmptyNode(3);
    PhaseSimulationMesh Box(EmptyNode);
    // Box.findNode(23,32,0).updateNode(12345);
    // Box.findNode(23,32,0).showNode();
    // Box.showProperties();
    EmptyNode.showNode();
=======
int main(){
    PhaseSimulationMesh Box(std::vector<long> {16,3,2});
    vector<PhaseNode> nodes(3);
    SpinodalSolver Solver;
    PhaseDomain NewDomain(nodes,Box);
    NewDomain.showBaseMeshProp();
    NewDomain.showGrainProp();
    Solver.showBaseMeshProp();
    Solver.showGrainProp();
>>>>>>> 90a66cfd1410eb594c1bec8af8f00122f719aabe
    return 0;
}