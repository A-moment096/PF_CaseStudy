#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

int main()
{
    PhaseNode EmptyNode;
    PhaseSimulationMesh Box(EmptyNode);
    Box.findNode(3,1,0).updateNode(12345);
    Box.findNode(3,1,0).showNode();
    Box.showNodeProp(1);
    return 0;
}