#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>

#include "SimulationMesh.hh"
#include "write_vtk_grid_values.cpp"

void write_vtk_grid_values(int nx,int ny,double dx,double dy,int istep,std::vector<std::vector<double>> data);