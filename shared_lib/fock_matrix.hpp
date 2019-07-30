#include <vector>
#include <Eigen/Dense>

# pragma once

using std::vector; using std::string;

typedef Eigen::MatrixXd matrix;

int atom(int ao_index, int orbitals_per_atom);

int orb_index(int ao_index, int orbitals_per_atom);

matrix calculate_fock_matrix(matrix hamiltonian_matrix, matrix interaction_matrix, matrix density_matrix, int orbitals_per_atom, double dipole);
