#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>
#include <utility>

#include "SimulationMesh.hh"

template<typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T> vec )
{
    for(auto ele : vec)out<<ele<<"";
    return out;
};

auto RunTimeCounter(std::chrono::_V2::system_clock::time_point _referTime){
    std::chrono::high_resolution_clock::time_point thisTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(thisTime-_referTime);
    std::cout<<"Time taken: "<<( double )duration.count()/1e6<<" seconds"<<std::endl;
    return ( double )duration.count()/1e6;
}
