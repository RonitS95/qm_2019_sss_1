class Name_Class: # Name of the class, and if there is nothing brackets, its an independent class
    def __init__(self, weight): # __init__ a constructor, self points to the object which calls the class, weight is just 
                                # an attribute of the class
        if isinstance(weight, float): # checks if the weight is a float or nor, if its not raises a TypeError
            self.weight = weight    
        else:
            raise TypeError("Weight needs to be a float")
    def __str__(self):  # If we dont have this __str__, then the class will return an address
        return ('Atomic weight is:', self.weight) # We refer to the attributes as self.attributes to that they can be 
                                                  # accessed outside of the class. 

Hydrogen = Name_Class(weight=1.008)               # Hydrogen: The self now refers to Hydrogen, and its assigns the weight to H 

#print(Hydrogen.weight)

Oxygen = Name_Class(weight=16.04)                 # self now refers to oxygen and assigns the weight to oxygen 


def molecular_weight(a,b):                       # A test function which takes the class attributes as arguments
    return a+b

def test_function():                             # A test function which has class attributes in its definition, with no input
    return  Hydrogen.weight


print(molecular_weight(2.0*Hydrogen.weight,Oxygen.weight)) # Both functions work
print(test_function())

#water = Name_Class(weight='18') # Raises an error

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#

class Forms_Molecule: # The usual class name
    def __init__(self, mol1, mol2): # class has two inputs, mol1 and mol2 and and an inside attribute, forms_bonds (empty str) 
        if isinstance(mol1, str):   # check for mol1 and mol2 to be str else raise an error 
            self.mol1 = mol1
        else: 
            raise TypeError("Molecule name needs to be a string")

        if isinstance(mol2, str):
            self.mol2 = mol2
        else:
            raise TypeError("Molecule name needs to be a string")

        self.forms_bond = ''  # declaration of forms_bond attribute

    def forms_bonds(self, bond_length): # definition of the inside attribute
        if isinstance(bond_length, float): # a function inside the class, self again points to the object which calls the 
            if bond_length > 1.5:          # function
                self.forms_bond = 'No'    # after it is defined, it needs to be assigned a value by the object outside.
            else:
                self.forms_bond = 'Yes'
        else:
            raise TypeError("Bond length needs to be a float")

    def __str__(self): # The usual str constructor
        return 'Molecule 1:' + self.mol1 + 'Molecule 2:' + self.mol2 + 'form bond? ' + self.forms_bond

Water = Forms_Molecule('Oxygen', 'Hydrogen')
Water.forms_bonds(1.3) # This is how the object assigns the inside attribute and below is how its accessed. 
print(Water.forms_bond)
