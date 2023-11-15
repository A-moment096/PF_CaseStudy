#include <fstream>
#include <vector>

void write_vtk_grid_values(int nx,int ny,double dx,double dy,int istep,std::vector<double> data)
{
	int nz = 1;
	int npoin = nx*ny*nz;

	char filename[128];
	sprintf(filename,"L:\\Programme\\C++\\PhaseFieldModelling\\CaseStudy_3\\output\\Result\\time_%06d.vtk", istep);

	std::ofstream outfile;
	outfile.open(filename);
	outfile<<"# vtk DataFile Version 2.0\n";
	outfile<<"time_10.vtk\n";
	outfile<<"ASCII\n";
	outfile<<"DATASET STRUCTURED_GRID\n";

	outfile<<"DIMENSIONS "<<nx<<"  "<<ny<<"  "<<nz<<"\n";
	outfile<<"POINTS "<<npoin<<"   float\n";
	double dumx,dumy,dumz;
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++)
		{
			dumx = i*dx;dumy = j*dy;dumz = 0.0;

			outfile<<dumx<<"   "<<dumy<<"   "<<dumz<<"\n";
		}

	outfile<<"POINT_DATA "<<npoin<<"\n";
	outfile<<"SCALARS CON  float  1\n";
	outfile<<"LOOKUP_TABLE default\n";

	for(int i=0;i<nx*ny;i++)
		outfile<<data.at(i)<<"\n";

	outfile.close();
	
}
