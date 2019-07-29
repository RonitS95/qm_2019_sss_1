#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "fock_matrix.hpp"

PYBIND11_MODULE(hartree_fock, m) // name of the output module, variable
{
    m.doc() = "This C++ library computes the Fock matrix in a Python environment.";
    m.def("atom", atom, "Gives back the atom number");
    m.def("orb_index", orb_index, "Gives back the orbital index");
    m.def("calculate_fock_matrix", calculate_fock_matrix, "Calculates the Fock matrix");
}
