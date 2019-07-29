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

int orb_index(int ao_index, int orbitals_per_atom)
{
    return ao_index % orbitals_per_atom;
}

int ao_index(int atom_p, string orb_p, vector<string> orbitals, int orbitals_per_atom)
{
    int p = atom_p * orbitals_per_atom;
    vector<string>::iterator itr = std::find(orbitals.begin(), orbitals.end(), orb_p);
    int dist;
    dist = std::distance(orbitals.begin(), itr);
    p = p + dist;

    return p;
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

float chi_on_atom(int o1, int o2, int o3, double dipole)
{
    if (o1 == o2 && o3 == 0)
        return 1.0;
    else if (o1 == o3 && (o3 > 0 && o3 <= 3) && o2 == 0)
        return dipole;
    else if (o2 == o3 && (o3 > 0 && o3 <= 3) && o1 == 0)
        return dipole;
    return 0.0;
}

matrix calculate_fock_matrix(matrix hamiltonian_matrix, matrix interaction_matrix, matrix density_matrix, vector<string> orbitals, int orbitals_per_atom, double dipole){
    int ndof = hamiltonian_matrix.row(0).size();
    matrix fock_matrix = hamiltonian_matrix;

    for (int p = 0; p < ndof; p++)
    {
        int at_p = atom(p, orbitals_per_atom);
        int orb_p = orb_index(p, orbitals_per_atom);
        for (int orb_q = 0; orb_q < orbitals_per_atom; orb_q++)
        {
            int q = orb_q + at_p * orbitals_per_atom;
            for (int orb_t = 0; orb_t < orbitals_per_atom; orb_t++)
            {
                int t = orb_t + at_p * orbitals_per_atom;
                float chi_pqt = chi_on_atom(orb_p, orb_q, orb_t, dipole);
                for (int r = 0; r < ndof; r++)
                {
                    int at_r = atom(r, orbitals_per_atom);
                    int orb_r = orb_index(r, orbitals_per_atom);
                    for (int orb_s = 0; orb_s < orbitals_per_atom; orb_s++)
                    {
                        int s = orb_s + at_r * orbitals_per_atom;
                        for (int orb_u = 0; orb_u < orbitals_per_atom; orb_u++)
                        {
                            int u = orb_u + at_r * orbitals_per_atom;
                            float chi_rsu = chi_on_atom(orb_r, orb_s, orb_u, dipole); 
                            fock_matrix(p,q) += 2.0 * chi_pqt * chi_rsu * interaction_matrix(t, u) * density_matrix(r, s);
                        }
                    }
                }
            }
        }
    }
    for (int p = 0; p < ndof; p++)
    {
        int at_p = atom(p, orbitals_per_atom);
        int orb_p = orb_index(p, orbitals_per_atom);
        for (int orb_s = 0; orb_s < orbitals_per_atom; orb_s++)
        {
            int s = orb_s + at_p * orbitals_per_atom;
            for (int orb_u = 0; orb_u < orbitals_per_atom; orb_u++)
            {
                int u = orb_u + at_p * orbitals_per_atom;
                float chi_psu = chi_on_atom(orb_p, orb_s, orb_u, dipole); 
                for (int q = 0; q < ndof; q++)
                {
                    int at_q = atom(q, orbitals_per_atom);
                    int orb_q = orb_index(q, orbitals_per_atom);
                    for (int orb_r = 0; orb_r < orbitals_per_atom; orb_r++)
                    {
                        int r = orb_r + at_q * orbitals_per_atom;
                        for (int orb_t = 0; orb_t < orbitals_per_atom; orb_t++)
                        {
                            int t = orb_t + at_q * orbitals_per_atom;
                            float chi_rqt = chi_on_atom(orb_r, orb_q, orb_t, dipole);
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
    matrix fock_matrix = calculate_fock_matrix(ham_test, int_test, dens_test, orbitals, orbitals_per_atom, dipole);

    cout << fock_matrix << endl;

    return 0;
}
