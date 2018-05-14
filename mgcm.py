import numpy as np  # knihovna NumPy
import sys  # pro ukončení práce při chybném souboru


class Calculate:
    def __init__(self, molecules, periodic_table):
        self.molecules, self.output, atom_parameter, not_errors_molecules = molecules, [], {}, 0
        for molecule in self.molecules:
            try:
                """příprava matic"""
                data_from_atoms, name, data_bond = [], molecule.name, []
                count = molecule.count_atoms
                degree_matrix = molecule.count_bond_matrix
                connectivity_matrix = molecule.bond_matrix
                identity_matrix = np.eye(count)
                """výpočet S = D - A + I"""
                simplified_matrix = degree_matrix - connectivity_matrix + identity_matrix
                pt_electronegativity = np.zeros((count, 1))
                try:
                    """příprava vektoru s hodnotami elektronegativity z periodické tabulky """
                    for i, atom in enumerate(molecule.atoms):
                        pt_electronegativity[i][0] = periodic_table[atom.element_symbol]
                        data_from_atoms.append((atom.element_symbol, atom.bond))
                except IndexError:
                    print("Can not calculate with ", name)
                    continue
                not_errors_molecules += 1  # součet spočítaných molekul
                """výpočet X = S^{-1} * X^{0}"""
                nk_electronegativity = np.linalg.solve(simplified_matrix, pt_electronegativity)
                """výpočet X - X^{0}"""
                deviation_away_pt = nk_electronegativity-pt_electronegativity
                multiple = 1
                """výpočet geometrického průměru"""
                for electroneg in pt_electronegativity:
                    multiple = multiple * electroneg
                geometric_mean = multiple**(1/count)
                """výpočet náboje"""
                charges = deviation_away_pt * (1/geometric_mean)
                self.output.append((name, count, data_from_atoms, charges))
            except KeyError:
                print("Something wrong with calculate")
                sys.exit()
        print("Program calculated {} molecules.".format(not_errors_molecules))

    def save_charges(self, file):
        data = self.output
        new_file = "result/" + file
        with open(new_file, "w") as f:
            for name, count, atoms, charges in data:
                print("{}\n{}".format(name, int(count)), file=f)
                for i, (atom, bond) in enumerate(atoms):
                    print("{0:6d}  {1:>2}{2} {3: f}".format(i + 1, atom, bond, float(charges[i])), file=f)
        print("Now you can find charge for each element in file {}".format(new_file))

    def give_result(self):
        return self.output
