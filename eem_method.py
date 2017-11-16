import sys
import argparse
import math
import numpy as np

from collections import Counter


class Atom:
    def __init__(self, number, element_symbol, bond, coordinate):
        self.element_symbol = element_symbol
        self.coordinate = coordinate
        self.bond = bond
        self.number = number

    def __str__(self):
        return str("{} {} {} {}".format(self.number, self.element_symbol, self.bond, self.coordinate))


class Molecule:

    def __init__(self, name, atoms, elements_count):
        self.name = name
        self.atoms = atoms
        self.elements_count = elements_count

    def __str__(self):
        return str("Molecule name: {}\n{}".format(self.name, self.atoms))


class MoleculesSet:
    def load_from_sdf(self, filename):
        self.filename = filename
        molecules = []
        try:
            with open(filename, "r") as fh:
                while True:
                    bond_data, atoms, elements, coordinate, line, elements_count = {}, [], [], [], fh.readline(),\
                                                                                   Counter()
                    if "" == line[0:1]:
                        self.molecules = molecules
                        return print("Load molecules from {}".format(filename))
                    name = (line[:].strip())
                    for i in range(2):
                            fh.readline()
                    line = fh.readline()
                    count_atoms, count_bonds = int(line[0:3]), int(line[3:6])
                    for i in range(1, count_atoms + 1):
                        line = fh.readline()
                        bond_data[i] = 0
                        elements.append(line[31:32])
                        coordinate.append((float(line[2:10]), float(line[12:20]), float(line[22:30])))
                    for i in range(count_bonds):
                        line = fh.readline()
                        first_atom, second_atom, bond = int(line[1:3]), int(line[3:6]), int(line[8:9])
                        bond_data[first_atom] = max(bond_data[first_atom], bond)
                        bond_data[second_atom] = max(bond_data[second_atom], bond)
                    for number in bond_data:
                        elements_count[(elements[number - 1], bond_data[number])] += 1
                        atoms.append((Atom(number, elements[number - 1], bond_data[number], coordinate[number - 1])))
                    while True:
                        line = fh.readline()
                        if "$$$$" in line:
                            molecules.append((Molecule(name, atoms, elements_count)))
                            break
        except IOError:
            print("Wrong file for molecules set! Try another file than {}".format(filename))
            sys.exit(1)

    def load_parameters(self, filename):
        try:
            with open(filename, "r") as fp:
                parameters = []
                for line in fp:
                    if " Kappa=" in line:
                        self.kappa = float(line[line.index("Kappa") + len("Kappa=\""):-3])
                    if "<Element" in line:
                        element = line[line.index("Name") + len("Name=\"")]
                        element_parameters = []
                    if "<Bond" in line:
                        parameter_a = int(line.index("A=") + len("A=\""))
                        parameter_b = int(line.index("B=") + len("B=\""))
                        element_parameters.append((float(line[parameter_a:parameter_a + 6]),
                                                  float(line[parameter_b:parameter_b + 6])))
                    if "</Element>" in line:
                        parameters.append((element, element_parameters))
                self.parameters = parameters
            return print("Load parameters for elements from {}".format(filename))
        except IOError:
            print("Wrong file for parameters! Try another file than").format(filename)
            sys.exit(2)

    def __str__(self):
        return str("{}".format(self.molecules))

    def calculate_atom_charge(self):
        atom_parameter, output = {}, []
        try:
            for atom, parameter in self.parameters:
                for bond, item in enumerate(parameter, start=1):
                    a_parameter, b_parameter = item
                    atom_parameter[atom, "A", bond] = a_parameter
                    atom_parameter[atom, "B", bond] = b_parameter
            for molecule in self.molecules:
                elements, bonds, numbers_elements, row_delete = [], [], [], []
                name = molecule.name
                for atom in molecule.atoms:
                    elements.append(atom.element_symbol)
                    bonds.append(atom.bond)
                    numbers_elements.append(atom.number)
                count = len(elements)
                distance = np.zeros((count + 1, count + 1))
                parameters_a = np.zeros((count + 1))
                for i, item in enumerate(elements):
                    try:
                        parameters_a[i] = - atom_parameter[item, "A", bonds[i]]
                    except KeyError:
                        parameters_a[i] = 0
                        continue
                for i, item in enumerate(elements):
                    try:
                        for j in range(i, count):
                                if i == j:
                                    distance[i, j] = atom_parameter[item, "B", bonds[i]]
                                else:
                                    distance[i, j] = distance[j, i] = get_distance(self.kappa,
                                                                                   molecule.atoms[i].coordinate,
                                                                                   molecule.atoms[j].coordinate)
                    except KeyError:
                        print("Missing parameters for {}. element ({}) in {}. Program did not count with this element."
                              "".format(i + 1, item, name))
                        row_delete.append(i)
                        continue
                distance[count, :] = 1
                distance[:, count] = -1
                distance[count, count] = 0
                if len(row_delete) > 0:
                    distance = np.delete(distance, row_delete, 1)
                    distance = np.delete(distance, row_delete, 0)
                    parameters_a = np.delete(parameters_a, row_delete, 0)
                    elements = [i for j, i in enumerate(elements) if j not in row_delete]
                    numbers_elements = [i for j, i in enumerate(numbers_elements) if j not in row_delete]
                charges = np.linalg.solve(distance, parameters_a)
                output.append((name, count, numbers_elements, elements, charges))
        except KeyError:
            print("Something wrong with calculate")
        output_to_file(self.filename, output)

    def get_statistic_from_set(self):
        elements, count_element = [], Counter()
        for molecule in self.molecules:
            for atom in molecule.atoms:
                elements.append(atom.element_symbol)
            count_element += molecule.elements_count
            print("Molecule:{}{}".format(molecule.name, "".join("\n{} = {}".format(atom, count) for atom,
                                                                count in molecule.elements_count.items())))
        elements_numbers = Counter(elements)
        print("Number of elements in whole set {}:".format(self.filename))
        for key in elements_numbers:
            print("{} = {}".format(key, elements_numbers[key]))
        print("Number of elements in whole set {} by bond:".format(self.filename))
        for key in sorted(count_element):
            print("{} = {}".format(key, count_element[key]))


def output_to_file(filename, data):
    new_file = "data_from_" + filename
    with open(new_file, "w") as f:
        for name, count, numbers, elements, charges in data:
            print("{}\n{}".format(name, count), file=f)
            for i, item in enumerate(elements):
                print("{0:6d}  {1:3} {2: f}".format(numbers[i], item, charges[i]), file=f)
            print("", file=f)
    print("Now you can find charge for each element in file {}".format(new_file))


def get_distance(kappa, first_coordinate, second_coordinate):
    x1, y1, z1 = first_coordinate
    x2, y2, z2 = second_coordinate
    return kappa / (math.sqrt(((x1 - x2) ** 2) + ((y1 - y2) ** 2) + ((z1 - z2) ** 2)))


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("SDF_file", help="Give a file with molecules (.sdf) ", type=str)
    parser.add_argument("Parameters_file", help="Give a file with parameters in xml format", type=str)
    args = parser.parse_args()
    mset = MoleculesSet()
    mset.load_from_sdf(args.SDF_file)
    mset.get_statistic_from_set()
    mset.load_parameters(args.Parameters_file)
    mset.calculate_atom_charge()


if __name__ == "__main__":
    main()
