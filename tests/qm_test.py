import numpy as np


class Nobel_Gas_model:
    def __init__(self, gas):
        if isinstance(gas, str): 
            self.gas = gas
        else:
            raise TypeError('Gas Name should be str')

        self.model_parameters = {}
        self.ionic_charge = 6 
        self.orbital_types = ['s', 'px', 'py', 'pz']
        self.p_orbitals = self.orbital_types[1:]
        self.orbitals_per_atom = len(self.orbital_types)
        self.vec = {'px': [1.0,0.0,0.0], 'py' : [0.0, 1.0, 0.0], 'pz' : [0.0, 0.0, 1.0]}
        self.orbital_occupations = {'s' : 0, 'px' : 1, 'py' : 1, 'pz' : 1 }
       
        if self.gas == 'Neon':
            self.model_parameters = {
                                     'coulomb_p': -0.010255409806855187,
                                     'coulomb_s': 0.4536486561938202,
                                     'dipole': 1.6692376991516769,
                                     'energy_p': -3.1186533988406335,
                                     'energy_s': 11.334912902362603,
                                     'r_hop': 2.739689713337267,
                                     'r_pseudo': 1.1800779720963734,
                                     't_pp1': -0.029546671673199854,
                                     't_pp2': -0.0041958662271044875,
                                     't_sp': 0.000450562836426027,
                                     't_ss': 0.0289251941290921,
                                     'v_pseudo': -0.015945813280635074
                                    }
        elif self.gas == 'Argon':
            self.model_parameters = { 
                                     'r_hop' : 3.1810226927827516, # hopping length scale
                                     't_ss' : 0.03365982238611262, # s-s hopping energy scale
                                     't_sp' : -0.029154833035109226, # s-p hopping energy scale
                                     't_pp1' : -0.0804163845390335, # 1st p-p hopping energy scale
                                     't_pp2' : -0.01393611496959445, # 2nd p-p hopping energy scale
                                     'r_pseudo' : 2.60342991362958, # pseudopotential length scale
                                     'v_pseudo' : 0.022972992186364977, # pseudopotential energy scale
                                     'dipole' : 2.781629275106456, # dipole strength of s-p transition
                                     'energy_s' : 3.1659446174413004, # onsite energy of s orbital
                                     'energy_p' : -2.3926873325346554, # onsite energy of p orbital
                                     'coulomb_s' : 0.3603533286088998, # Coulomb self-energy of monopole
                                     'coulomb_p' : -0.003267991835806299  # Coulomb self-energy of dipole
                                    }
        else:
            return 'Enter a valid gas'

        #self.orb = self.orb(ao_index)
       # self.atom = self.atom(ao_index)
       # self.ao_index = self.ao_index(atom_p, orb_p)

    def orb(self, ao_index):
        orb_index = ao_index % self.orbitals_per_atom
        return self.orbital_types[orb_index]
    def atom(self, ao_index):
        return ao_index // self.orbitals_per_atom
    def ao_index(atom_p, orb_p):
        p = atom_p * orbitals_per_atom
        p += self.orbital_types.index(orb_p)
        return p
    
    
    def __str__(self):
        return 'Model Parameters for' + self.gas 


# Need to define number of atoms here

# Need to define coordinates here since HF class initializes with those! 

class Hartree_Fock:
    def __init__(self, Gas_Model, atomic_coordinates):
        self.Gas_Model = Nobel_Gas_model(Gas_Model)
        self.atomic_coordinates = atomic_coordinates
        self.ndof = len(self.atomic_coordinates)*self.Gas_Model.orbitals_per_atom
        



    


    def hopping_energy(self, orb1, orb2, r12, model_parameters):
        r12_rescaled = r12 / model_parameters['r_hop']
        r12_length = np.linalg.norm(r12_rescaled)
        ans = np.exp(1.0 - r12_length**2)
        if orb1 == 's' and orb2 == 's':
            ans *= model_parameters['t_ss']
        if orb1 == 's' and orb2 in self.Gas_Model.p_orbitals:
            ans *= np.dot(self.Gas_Model.vec[orb2], r12_rescaled)*model_parameters['t_sp']
        if orb1 == 'p' and orb1 in self.Gas_Model.p_orbitals:
            ans *= - np.dot(self.Gas_Model.vec[orb1], r12_rescaled)*model_parameters['t_sp']
        if orb1 in self.Gas_Model.p_orbitals and orb2 in self.Gas_Model.p_orbitals:
            ans *= ( (r12_length**2) * np.dot(self.Gas_Model.vec[orb1], self.Gas_Model.vec[orb2]) * model_parameters['t_pp2']
                 - np.dot(self.Gas_Model.vec[orb1], r12_rescaled) * np.dot(self.Gas_Model.vec[orb2], r12_rescaled)
                 * ( model_parameters['t_pp1'] + model_parameters['t_pp2'] ) )
        return ans


    def __str__(self):
        return 'ndof:' + str(self.ndof)
# Testing-------------

system = Nobel_Gas_model('Argon')

print(system.model_parameters['r_hop'])

for index in range(2*system.orbitals_per_atom):
    print('index: ', index, 'atom', system.atom(index), 'orbital', system.orb(index))

coord = [[1,0,0], [5,0,0]]
calc = Hartree_Fock('Argon', coord)
print(calc.__dict__)

print('hopping tests:')
hop_vector = np.array([system.model_parameters['r_hop'],0,0])
print('s s =', calc.hopping_energy('s','s',hop_vector,system.model_parameters), system.model_parameters['t_ss'])
print('s px =', calc.hopping_energy('s','px',hop_vector,system.model_parameters), system.model_parameters['t_sp'])
print('px px =', calc.hopping_energy('px','px',hop_vector,system.model_parameters), system.model_parameters['t_pp1'])
print('s py =', calc.hopping_energy('s','py',hop_vector,system.model_parameters))
print('py px =', calc.hopping_energy('py','px',hop_vector,system.model_parameters))
print('py py =', calc.hopping_energy('py','py',hop_vector,system.model_parameters), system.model_parameters['t_pp2'])
# Testing--------------
