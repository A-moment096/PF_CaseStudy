#pragma once
#ifndef PFMTOOLS_HH
#define PFMTOOLS_HH

#include <iostream>
#include <chrono>
#include <vector>

namespace PFMTools
{

    /*************************************************************/

    double threshold(double &val, double min, double max){
        return val>max ? val = max : (val<min ? val = min : val = val);
    }
    
    double threshold(double &&val, double min, double max){
        return val>max ? val = max : (val<min ? val = min : val = val);
    }

    double threshold(double &val){
        return val>0.999999 ? val = 0.999999 : (val<0.000001 ? val = 0.000001 : val = val);
    }

    double threshold(double &&val){
        return val>0.999999 ? val = 0.999999 : (val<0.000001 ? val = 0.000001 : val = val);
    }

    /*************************************************************/

    template<typename T>
    std::ostream& operator<<(std::ostream& out, const std::vector<T> vec )
    {
        for(auto ele : vec)out<<ele<<"";
        return out;
    };
    
    /*************************************************************/

    using TimePoint = std::chrono::_V2::system_clock::time_point;
    
    TimePoint now(){
        return std::chrono::high_resolution_clock::now();
    }

    double RunTimeCounter(std::chrono::_V2::system_clock::time_point _referTime, int identifier){
        std::chrono::high_resolution_clock::time_point thisTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(thisTime-_referTime);
        // std::cout<<"\tTime taken: "<<( double )duration.count()/1e6<<" seconds"<<std::endl;
        return ( double )duration.count()/1e6;
    }

    void RunTimeCounter(std::chrono::_V2::system_clock::time_point _referTime){
        std::chrono::high_resolution_clock::time_point thisTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(thisTime-_referTime);
        std::cout<<"\tTime taken: "<<( double )duration.count()/1e6<<" seconds"<<std::endl;
        return ;
    }
} // namespace PFMTools



#endif