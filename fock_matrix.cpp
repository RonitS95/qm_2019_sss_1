#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <algorithm>

using std::cout; using std::endl; using std::vector; using std::string;

typedef Eigen::MatrixXd matrix;

int atom(int ao_index, int orbitals_per_atom)
{
    return ao_index / orbitals_per_atom;
}

int ao_index(int atom_p, string orb_p, vector<string> orbitals, int orbitals_per_atom)
{
    int p = atom_p * orbitals_per_atom;
    //std::vector<int>::iterator
    auto itr = std::find(orbitals.begin(), orbitals.end(), orb_p);
    int dist;
    dist = std::distance(orbitals.begin(), itr);
    p = p + dist;

    return p;
}

string orb(int ao_index, vector<string> orbitals, int orbitals_per_atom)
{
    int orb_index = ao_index % orbitals_per_atom;
    
    return orbitals.at(orb_index);
}

template <typename T, typename U>
bool in_vector(U elem, const vector<T>& vec)
{
    return (std::find(vec.begin(), vec.end(), elem) != vec.end()) ? true : false;
}

template <typename T>
vector<T> slice(vector<T>& vec, int i, int f)
{
    vector<T> sliced(f - i + 1);
    std::copy(vec.begin() + i, vec.begin() + f + 1, sliced.begin());
    return sliced;
}

float chi_on_atom(string o1, string o2, string o3, vector<string> orbitals, double dipole)
{
    if (o1 == o2 && o3 == "s")
        return 1.0;
    else if (o1 == o3 && in_vector(o3, slice(orbitals, 1, orbitals.size() - 1)) && o2 == "s")
        return dipole;
    else if (o2 == o3 && in_vector(o3, slice(orbitals, 1, orbitals.size() - 1)) && o1 == "s")
        return dipole;
    return 0.0;
}

matrix calculate_fock_matrix(matrix hamiltonian_matrix, matrix interaction_matrix, matrix density_matrix, vector<string> orbitals, int orbitals_per_atom, double dipole)
{
    int ndof = hamiltonian_matrix.row(0).size();
    matrix fock_matrix = hamiltonian_matrix;

    for (int p = 0; p < ndof; p++)
    {
        for (int orb_q = 0; orb_q < orbitals_per_atom; orb_q++)
        {
            int q = ao_index(atom(p, orbitals_per_atom), orbitals.at(orb_q), orbitals, orbitals_per_atom);
            for (int orb_t = 0; orb_t < orbitals_per_atom; orb_t++)
            {
                int t = ao_index(atom(p, orbitals_per_atom), orbitals.at(orb_t), orbitals, orbitals_per_atom);
                float chi_pqt = chi_on_atom(orb(p, orbitals, orbitals_per_atom), orbitals.at(orb_q), orbitals.at(orb_t), orbitals, dipole);
                for (int r = 0; r < ndof; r++)
                {
                    for (int orb_s = 0; orb_s < orbitals_per_atom; orb_s++)
                    {
                        int s = ao_index(atom(r, orbitals_per_atom), orbitals.at(orb_s), orbitals, orbitals_per_atom);
                        for (int orb_u = 0; orb_u < orbitals_per_atom; orb_u++)
                        {
                            int u = ao_index(atom(r, orbitals_per_atom), orbitals.at(orb_u), orbitals, orbitals_per_atom);
                            float chi_rsu = chi_on_atom(orb(r, orbitals, orbitals_per_atom), orbitals.at(orb_s), orbitals.at(orb_u), orbitals, dipole);
                            fock_matrix(p,q) += 2.0 * chi_pqt * chi_rsu * interaction_matrix(t, u) * density_matrix(r, s);
                        }
                    }
                }
            }
        }
    }
    for (int p = 0; p < ndof; p++)
    {
        for (int orb_s = 0; orb_s < orbitals_per_atom; orb_s++)
        {
            int s = ao_index(atom(p, orbitals_per_atom), orbitals.at(orb_s), orbitals, orbitals_per_atom);
            for (int orb_u = 0; orb_u < orbitals_per_atom; orb_u++)
            {
                int u = ao_index(atom(p, orbitals_per_atom), orbitals.at(orb_u), orbitals, orbitals_per_atom);
                float chi_psu = chi_on_atom(orb(p, orbitals, orbitals_per_atom), orbitals.at(orb_s), orbitals.at(orb_u), orbitals, dipole);
                for (int q = 0; q < ndof; q++)
                {
                    for (int orb_r = 0; orb_r < orbitals_per_atom; orb_r++)
                    {
                        int r = ao_index(atom(q, orbitals_per_atom), orbitals.at(orb_r), orbitals, orbitals_per_atom);
                        for (int orb_t = 0; orb_t < orbitals_per_atom; orb_t++)
                        {
                            int t = ao_index(atom(q, orbitals_per_atom), orbitals.at(orb_t), orbitals, orbitals_per_atom);
                            float chi_rqt = chi_on_atom(orb(r, orbitals, orbitals_per_atom), orb(q, orbitals, orbitals_per_atom), orbitals.at(orb_t), orbitals, dipole);
                            fock_matrix(p,q) -= chi_rqt * chi_psu * interaction_matrix(t, u) * density_matrix(r, s);
                        }
                    }
                }
            }
        }
    }
    return fock_matrix;
}

int main()
{
    Eigen::MatrixXd ham_test(8,8);
    ham_test << 2.31745554e+00, -1.41367040e-01, -1.88489387e-01, -2.35611734e-01,
    6.53808346e-04,  5.34076767e-04,  7.12102356e-04,  8.90127945e-04, -1.41367040e-01,
    -3.24117641e+00,  0.00000000e+00,  0.00000000e+00, -5.34076767e-04,  2.92479412e-04,
    2.17340049e-03,  2.71675062e-03,  -1.88489387e-01,  0.00000000e+00, -3.24117641e+00,
    0.00000000e+00,  -7.12102356e-04,  2.17340049e-03,  1.56029637e-03,  3.62233416e-03,
    -2.35611734e-01,  0.00000000e+00,  0.00000000e+00, -3.24117641e+00, -8.90127945e-04,
    2.71675062e-03,  3.62233416e-03,  3.19034674e-03,  6.53808346e-04, -5.34076767e-04,
    -7.12102356e-04, -8.90127945e-04,  2.31745554e+00,  1.41367040e-01,  1.88489387e-01,
    2.35611734e-01,  5.34076767e-04,  2.92479412e-04,  2.17340049e-03,  2.71675062e-03,
    1.41367040e-01, -3.24117641e+00,  0.00000000e+00,  0.00000000e+00,  7.12102356e-04,
    2.17340049e-03,  1.56029637e-03,  3.62233416e-03,  1.88489387e-01,  0.00000000e+00,
    -3.24117641e+00,  0.00000000e+00,  8.90127945e-04,  2.71675062e-03,  3.62233416e-03,
    3.19034674e-03,  2.35611734e-01,  0.00000000e+00,  0.00000000e+00, -3.24117641e+00;

    Eigen::MatrixXd int_test(8,8);
    int_test << 3.60353329e-01,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
   1.41421356e-01, -8.48528137e-03, -1.13137085e-02, -1.41421356e-02,
    0.00000000e+00, -3.26799184e-03,  0.00000000e+00,  0.00000000e+00,
   8.48528137e-03,  1.30107648e-03, -2.03646753e-03, -2.54558441e-03,
    0.00000000e+00,  0.00000000e+00, -3.26799184e-03,  0.00000000e+00,
   1.13137085e-02, -2.03646753e-03,  1.13137085e-04, -3.39411255e-03,
    0.00000000e+00,  0.00000000e+00,  0.00000000e+00, -3.26799184e-03,
   1.41421356e-02, -2.54558441e-03, -3.39411255e-03, -1.41421356e-03,
    1.41421356e-01,  8.48528137e-03,  1.13137085e-02,  1.41421356e-02,
   3.60353329e-01,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
   -8.48528137e-03,  1.30107648e-03, -2.03646753e-03, -2.54558441e-03,
   0.00000000e+00, -3.26799184e-03,  0.00000000e+00,  0.00000000e+00,
   -1.13137085e-02, -2.03646753e-03,  1.13137085e-04, -3.39411255e-03,
   0.00000000e+00,  0.00000000e+00, -3.26799184e-03,  0.00000000e+00,
   -1.41421356e-02, -2.54558441e-03, -3.39411255e-03, -1.41421356e-03,
   0.00000000e+00,  0.00000000e+00,  0.00000000e+00, -3.26799184e-03;

    Eigen::MatrixXd dens_test(8,8);
    dens_test << 0., 0., 0., 0., 0., 0., 0., 0.,
 0., 1., 0., 0., 0., 0., 0., 0.,
 0., 0., 1., 0., 0., 0., 0., 0.,
 0., 0., 0., 1., 0., 0., 0., 0.,
 0., 0., 0., 0., 0., 0., 0., 0.,
 0., 0., 0., 0., 0., 1., 0., 0.,
 0., 0., 0., 0., 0., 0., 1., 0.,
 0., 0., 0., 0., 0., 0., 0., 1.;

    vector<string> orbitals{"s", "px" ,"py", "pz"};
    int orbitals_per_atom = orbitals.size();
    double dipole(2.781629275106456);

    string test_str("px");

    int test_atom = atom(2, orbitals_per_atom);
    int test_ao_index = ao_index(2, "px", orbitals, orbitals_per_atom);
    string orb_test = orb(test_ao_index, orbitals, orbitals_per_atom);
    bool in_vector_test = in_vector("px", orbitals);
    vector<string> test_slice = slice(orbitals, 1, 2);
    float chi_test = chi_on_atom("px", "py", "px", orbitals, dipole);
    cout << orbitals_per_atom << " " << test_atom << " " << test_ao_index << " " << orb_test << " " << in_vector_test << " " << chi_test << endl;

     matrix fock_matrix2 = calculate_fock_matrix(ham_test, int_test, dens_test, orbitals, orbitals_per_atom, dipole);

     cout << fock_matrix2 << endl;

    return 0;
}
