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
        for(auto ele : vec)out<<ele<<" ";
        return out;
    };
    
    /*************************************************************/

    using TimePoint = std::chrono::_V2::system_clock::time_point;
    
    TimePoint now(){
        return std::chrono::high_resolution_clock::now();
    }

    double RunTimeCounter(std::chrono::_V2::system_clock::time_point _referTime, bool isprint){
        std::chrono::high_resolution_clock::time_point thisTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(thisTime-_referTime);
        if(isprint){
            std::cout<<"\tTime taken: "<<( double )duration.count()/1e6<<" seconds"<<std::endl;
            return 0;
        }
        return ( double )duration.count()/1e6;
    }

    /*************************************************************/
    //!!! 2D Methods Only !!!//
    std::vector<int> readCSV(std::string _path_name){
        std::vector<int> datas;
        std::ifstream csvfile;
        csvfile.open(_path_name);
        std::string line;
        while (std::getline(csvfile,line)){
            std::stringstream subline(line);
            while(std::getline(subline,line,',')){
                datas.push_back(std::stoi(line));
            }
        }
        csvfile.close();
        return datas;
    }

} // namespace PFMTools



#endif