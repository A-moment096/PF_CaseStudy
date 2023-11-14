#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

int main()
{
    PhaseNode EmptyNode(1);
    std::vector<int> Size{8,8,3};
    PhaseSimulationMesh Box({8,8,3},EmptyNode);
    for(int i = 0; i<Box.getNumNodes(); i++){
        Box.findNode(i).updateNode(WHICHPARA::CON,0,double(i));
    }
    Box.findNode(0).Backward->showNode();
    Box.findNode(0).Forward->showNode();
    Box.findNode(0).Up->showNode();
    Box.findNode(0).Down->showNode();
    Box.findNode(0).Left->showNode();
    Box.findNode(0).Right->showNode();
    
    
    return 0;
}