#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

int main(){
    PhaseNode EmptyNode;
    PhaseSimulationMesh Box(std::vector<long> {16,3,2},EmptyNode);
    Box.findNode(3,1,0).updateNode(std::vector<long double> {257,10.4,9});
    Box.findNode(3,1,0).showNode();
    Box.showNodesProp(1);
    return 0;
}