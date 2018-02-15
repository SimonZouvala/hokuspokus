import math
import numpy as np


class Calculate:
    def __init__(self, molecules, parameters):
        self.parameters, self.molecules, self.output, atom_parameter, error_molecules = parameters, molecules, [], {}, 0
        try:
            kappa, yes_type, parameters = self.parameters
            for element, type_bond, parameter in parameters:
                a_parameter, b_parameter = parameter[0]
                atom_parameter[element, "A", type_bond] = a_parameter
                atom_parameter[element, "B", type_bond] = b_parameter
            for molecule in self.molecules:
                data_from_atoms, name = [], molecule.name
                count = molecule.count_atoms
                distance = np.zeros((count + 1, count + 1))
                parameters_a = np.zeros((count + 1))
                try:
                    for i, atom in enumerate(molecule.atoms):
                        data_from_atoms.append((atom.element_symbol, atom.number, atom.bond))
                        if not yes_type:
                            atom.bond = 1
                        parameters_a[i] = - atom_parameter[atom.element_symbol, "A", atom.bond]
                        for j in range(i, count):
                            if i == j:
                                distance[i, j] = atom_parameter[atom.element_symbol, "B", atom.bond]
                            else:
                                distance[i, j] = distance[j, i] = get_distance(kappa,
                                                                               molecule.atoms[i].coordinate,
                                                                               molecule.atoms[j].coordinate)
                    error_molecules += 1
                except KeyError:
                    # print("Missing parameters for {}. element {}({}) in {}. Program did not count with this "
                    #     "element.".format(atom.number, atom.element_symbol, atom.bond, name))
                    self.output.append((name, "error", atom.element_symbol, 0, atom_parameter))
                    continue
                distance[count, :] = 1
                distance[:, count] = -1
                distance[count, count] = 0
                charges = np.linalg.solve(distance, parameters_a)
                self.output.append((name, count, data_from_atoms, charges, atom_parameter))
            print("Program calculated {} molecules.".format(error_molecules))
        except KeyError:
            print("Something wrong with calculate")

    def save_charges(self, new_file):
        data_for_graph, data = [], self.output
        self.new_file = "result/" + new_file
        with open(self.new_file, "w") as f:
            for name, count, atoms, charges, parameters in data:
                try:
                    print("{}\n{}".format(name, int(count)), file=f)
                    for i, atom in enumerate(atoms):
                        element, number, bond = atom
                        print("{0:6d}  {1:>2}{2} {3: f}".format(number, element, bond, charges[i]), file=f)
                    data_for_graph.append((name, count, atoms, charges, parameters))
                except ValueError:
                    continue
        print("Now you can find charge for each element in file {}".format(self.new_file))

    def give_result(self):
        return self.output


def get_distance(kappa, first_coordinate, second_coordinate):
    x1, y1, z1 = first_coordinate
    x2, y2, z2 = second_coordinate
    return kappa / (math.sqrt(((x1 - x2) ** 2) + ((y1 - y2) ** 2) + ((z1 - z2) ** 2)))
