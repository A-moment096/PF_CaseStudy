#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>

<<<<<<< HEAD
#include "PhaseSimulationMesh.hh"
=======
#include "SpinodalSolver.hh"
>>>>>>> 90a66cfd1410eb594c1bec8af8f00122f719aabe

void write_vtk_grid_values(int nx,int ny,double dx,double dy,int istep,std::vector<std::vector<double>> data);