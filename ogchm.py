import math
import numpy as np

from collections import Counter

class Calculate:
    def __init__(self, molecules, tb_el, tb_hard, covalent_radii):
        self.molecules, self.output, atom_parameter, error_molecules = molecules, [], {}, 0
        tb_electronegativity, tb_hardness = tb_el, tb_hard
        try:
            for molecule in self.molecules:
                data_from_atoms, name, data_bond = [], molecule.name, []
                count = molecule.count_atoms
                degree_matrix = molecule.count_bond_matrix
                connectivity_matrix = molecule.bond_matrix
                identity_matrix = np.eye(count)
                print("D = Degree matrix (OK) ")
                print(degree_matrix)
                print("A = Connectivity matrix (OK)")
                print(connectivity_matrix)
                print("I = Identity matrix (OK)")
                print(identity_matrix)
                simplified_matrix = degree_matrix - connectivity_matrix + identity_matrix
                print("S = Simplified matrix (OK)")
                print(simplified_matrix)
                print("X° = table electronegativity (OK)")
                print(tb_electronegativity)
                print("n = table Hardness")
                print(tb_hardness)
                nk_electronegativity = np.linalg.solve(simplified_matrix, tb_electronegativity)
                print("X = Unknown electronegativities")
                print(nk_electronegativity)
                deviation_away_pt = nk_electronegativity-tb_electronegativity
                print("X - X°")
                print(deviation_away_pt)
                denominator = 0.0
                numerator = 0.0
                for i in range(deviation_away_pt.shape[0]):
                    for j in range(deviation_away_pt.shape[0]):
                        denominator += deviation_away_pt[j]/tb_hardness[j]
                        numerator += (deviation_away_pt[j]**2)/((tb_hardness[j]**3) * covalent_radii[j])
                        print(covalent_radii[j])
                dm = denominator/numerator
                print(deviation_away_pt.shape[0])
                charges_electrons = (deviation_away_pt/tb_hardness)-(((deviation_away_pt**2)*dm)/(
                    (tb_hardness**3)*covalent_radii))
                shift = 0
                charge_elements = Counter()
                for element, number in molecule.elements_count:
                    for index in range(molecule.elements_count[element, number]):
                        charge_elements[element, number] += float(charges_electrons[index + shift])
                    shift += molecule.elements_count[element, number]
                    data_from_atoms.append((element, number))
                print(charge_elements)
                self.output.append((name, count, data_from_atoms, charge_elements))
        except KeyError:
            print("Something wrong with calculate")

    def save_charges(self, new_file):
        data = self.output
        print("255sssss-----")
        self.new_file = "result/" + new_file
        with open(self.new_file, "w") as f:
            for name, count, atoms, charges in data:
                print(name)
                print("{}\n{}".format(name, int(count)), file=f)
                for index, (atom, bond) in enumerate(atoms, 1):
                    print(atom, bond, charges[atom, index])
                    print("{0:<1} {1:<3} {2: 5f}".format(atom, bond, float(charges[atom, index])), file=f)
        print("Now you can find charge for each element in file {}".format(self.new_file))

    def give_result(self):
        return self.output


def get_distance(kappa, first_coordinate, second_coordinate):
    x1, y1, z1 = first_coordinate
    x2, y2, z2 = second_coordinate
    return kappa / (math.sqrt(((x1 - x2) ** 2) + ((y1 - y2) ** 2) + ((z1 - z2) ** 2)))

