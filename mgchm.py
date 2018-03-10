import numpy as np


class Calculate:
    def __init__(self, molecules, periodic_table):
        self.molecules, self.output, atom_parameter, error_molecules = molecules, [], {}, 0
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
                geometric_mean = np.sum(nk_electronegativity)/count
                charges = deviation_away_pt * (1/geometric_mean)
                print("Q = Charges")
                print(charges)
                self.output.append((name, count, data_from_atoms, charges))
        except KeyError:
            print("Something wrong with calculate")

    def save_charges(self, new_file):
        data = self.output
        self.new_file = "result/" + new_file
        with open(self.new_file, "w") as f:
            for name, count, atoms, charges in data:
                print("{}\n{}".format(name, int(count)), file=f)
                for i, (atom, bond) in enumerate(atoms):
                    print(atom, bond, charges[i])
                    print("{0:<1} {1:<3} {2: 5f}".format(atom, bond, float(charges[i])), file=f)
        print("Now you can find charge for each element in file {}".format(self.new_file))

    def give_result(self):
        return self.output
