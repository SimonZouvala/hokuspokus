import numpy as np
import math


class Calculate:
    def __init__(self, molecules, periodic_table):
        self.molecules, self.output, atom_parameter = molecules, [], {}
        try:
            for molecule in self.molecules:
                data_from_atoms, name, data_bond = [], molecule.name, []
                count = molecule.count_atoms
                degree_matrix = molecule.count_bond_matrix
                connectivity_matrix = molecule.bond_matrix
                identity_matrix = np.eye(count)
                simplified_matrix = degree_matrix - connectivity_matrix + identity_matrix
                pt_electronegativity = np.zeros((count, 1))
                for i, atom in enumerate(molecule.atoms):
                    pt_electronegativity[i][0] = periodic_table[atom.element_symbol]
                    data_from_atoms.append((atom.element_symbol, atom.bond))
                nk_electronegativity = np.linalg.solve(simplified_matrix, pt_electronegativity)
                deviation_away_pt = nk_electronegativity-pt_electronegativity
                multiple = 1
                for electroneg in pt_electronegativity:
                    multiple = multiple * electroneg
                geometric_mean = multiple**(1/count)
                charges = deviation_away_pt * (1/geometric_mean)
                self.output.append((name, count, data_from_atoms, charges))
                print("Program calculated {} molecules.".format(len(molecules)))
        except KeyError:
            print("Something wrong with calculate")

    def save_charges(self, file):
        data = self.output
        new_file = "result/" + file
        with open(new_file, "w") as f:
            for name, count, atoms, charges in data:
                print(name)
                print("{}\n{}".format(name, int(count)), file=f)
                for i, (atom, bond) in enumerate(atoms):
                    print(atom, bond, charges[i])
                    print("{0:6d}  {1:>2}{2} {3: f}".format(i + 1, atom, bond, float(charges[i])), file=f)
        print("Now you can find charge for each element in file {}".format(new_file))

    def give_result(self):
        return self.output
