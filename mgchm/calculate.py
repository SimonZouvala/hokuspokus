import math
import numpy as np
from numpy.linalg import inv


class Calculate:
    def __init__(self, molecules, periodic_table):
        self.molecules, self.output, atom_parameter, error_molecules = molecules, [], {}, 0
        try:
            for molecule in self.molecules:
                data_from_atoms, name = [], molecule.name
                count = molecule.count_atoms
                degree_matrix = molecule.count_bond_matrix
                connectivity_matrix = molecule.bond_matrix
                identity_matrix = np.eye(count)
                simplified_matrix = degree_matrix - connectivity_matrix + identity_matrix
                pt_electronegativity = np.zeros((count, 1))
                for i, atom in enumerate(molecule.atoms):
                    pt_electronegativity[i][0] = periodic_table[atom.element_symbol]
                    data_from_atoms.append(atom.element_symbol)
                nk_electronegativity = np.linalg.solve(inv(simplified_matrix), pt_electronegativity)
                geometric_mean = np.sum(nk_electronegativity)/count
                charges = nk_electronegativity * (1/geometric_mean)
                self.output.append((name, count, data_from_atoms, charges))
        except KeyError:
            print("Something wrong with calculate")

    def save_charges(self, new_file):
        data = self.output
        self.new_file = "result/" + new_file
        with open(self.new_file, "w") as f:
            for name, count, atoms, charges in data:
                print("{}\n{}".format(name, int(count)), file=f)
                for i, atom in enumerate(atoms):
                    print(atom)
                    print(charges[i])
                    print("{0:5} {1: 5f}".format(atom, float(charges[i])), file=f)
        print("Now you can find charge for each element in file {}".format(self.new_file))

    def give_result(self):
        return self.output


def get_distance(kappa, first_coordinate, second_coordinate):
    x1, y1, z1 = first_coordinate
    x2, y2, z2 = second_coordinate
    return kappa / (math.sqrt(((x1 - x2) ** 2) + ((y1 - y2) ** 2) + ((z1 - z2) ** 2)))
