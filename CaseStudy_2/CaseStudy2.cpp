#include "../include/PFM.hh"

const double mobility{5.0}, kappa{0.1}, A{1.0}, B{1.0};
const double dtime{0.005};

double dfdphi(double phi, double phi2_sum, double A, double B) {
    return -1.0 * A * phi + B * phi * phi * phi + 2 * (phi * phi2_sum - phi * phi * phi);
}

int main() {
    PFMTools::TimePoint start = PFMTools::now(), dur = PFMTools::now();
    std::string _path(toVTK_Path("./"));

    MeshNode node;
    node.addEnt(WHICHPARA::PHSFRAC, 1);
    SimulationMesh mesh({64, 64, 1}, {0.5, 0.5, 1}, dtime, std::move(node));

    std::cout << "ready" << std::endl;

    int nstep{5000}, nprint{100};
#pragma omp parallel for collapse(2)
    for (int i = 0; i < mesh.MeshX; i++) {
        for (int j = 0; j < mesh.MeshY; j++) {
            if ((i - mesh.MeshX / 2) * (i - mesh.MeshX / 2) + (j - mesh.MeshY / 2) * (j - mesh.MeshY / 2) <= 14.0 * 14.0) {
                mesh.updateNodeVal(WHICHPARA::PHSFRAC, {i, j, 0}, 0, 1.0);
                mesh.updateNodeVal(WHICHPARA::PHSFRAC, {i, j, 0}, 1, 0.0);
            } else {
                mesh.updateNodeVal(WHICHPARA::PHSFRAC, {i, j, 0}, 1, 1.0);
                mesh.updateNodeVal(WHICHPARA::PHSFRAC, {i, j, 0}, 0, 0.0);
            }
        }
    }

    std::cout << "ready to go" << std::endl;

    for (int istep{0}; istep <= nstep; istep++) {

        mesh.Laplacian(WHICHPARA::PHSFRAC);

#pragma omp parallel for
        for (auto &node : mesh.SimuNodes) {
            for (int phi_index{0}; phi_index < 2; phi_index++) {

                double bulk_contribution{dfdphi(node.Phs_Node.getVal(phi_index), node.sumPhsFrac2(), A, B)};
                double interphase_contribution{-1 * kappa * node.getLap(WHICHPARA::PHSFRAC, phi_index)};
                double dFdphi = -1.0 * mobility * (bulk_contribution + interphase_contribution);

                node.Phs_Node.updateDVal(phi_index, dFdphi);
            }
        }

        mesh.iterateVal(WHICHPARA::PHSFRAC);

        if (istep % nprint == 0) {
            mesh.outVTKFilehead(_path, istep);
            mesh.outVTKAll(_path, WHICHPARA::PHSFRAC, istep);
            std::cout << "Done Step: " << istep;
            PFMTools::RunTimeCounter(dur, true);
            dur = PFMTools::now();
        }
    }

    PFMTools::RunTimeCounter(start, true);
}
