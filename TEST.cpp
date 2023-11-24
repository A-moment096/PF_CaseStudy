#include "PFM.hh"

using std::cout;
using std::endl;
using std::vector;

int main(){
    auto start =  std::chrono::high_resolution_clock:: now();
    double c = 12;
    cout<< c*4*c*-2*c-1<<endl;
    cout<<c*4*c-2*c-1<<endl;
    auto stop =  std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start);
    cout<<"\nTime taken by programme: "<<(double) duration.count() / 1e6 << " seconds"<<endl;

    // system("pause");
    return 0;
}

