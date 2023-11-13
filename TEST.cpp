#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

int main()
{
    PhaseNode EmptyNode(3);
    PhaseSimulationMesh Box(EmptyNode);
    // Box.findNode(23,32,0).updateNode(12345);
    // Box.findNode(23,32,0).showNode();
    // Box.showProperties();
    EmptyNode.showNode();
    return 0;
}