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

    def __init__(self, name, count_atoms, atoms, elements_count):
        self.name = name
        self.count_atoms = count_atoms
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
                            molecules.append((Molecule(name, count_atoms,  atoms, elements_count)))
                            break
        except IOError:
            print("Wrong file for molecules set! Try another file than {}".format(filename))
            sys.exit(1)

    def load_parameters(self, filename):
        self.yes_type = False
        try:
            self.file_parameters = filename
            with open(filename, "r") as fp:
                parameters = []
                for line in fp:
                    if " Kappa=" in line:
                        self.kappa = float(line[line.index("Kappa") + len("Kappa=\""):-3])
                    if "<Element" in line:
                        element = line[line.index("Name") + len("Name=\""):-3].strip()
                    if "<Bond" in line:
                        element_parameters = []
                        if "Type" in line:
                            self.yes_type = True
                            type_bond = int(line[line.index("Type") + len("Type=\"")])
                        else:
                            type_bond = 1
                        parameter_a = int(line.index("A=") + len("A=\""))
                        parameter_b = int(line.index("B=") + len("B=\""))
                        element_parameters.append((float(line[parameter_a:parameter_a + 6]),
                                                  float(line[parameter_b:parameter_b + 6])))
                        parameters.append((element, type_bond, element_parameters))
                self.parameters = parameters
            return print("Load parameters for elements from {}".format(filename))
        except IOError:
            print("Wrong file for parameters! Try another file than").format(filename)
            sys.exit(2)
        get_statistic_from_parameters

    def calculate_atom_charge(self):
        atom_parameter, output = {}, []
        try:
            for element, type_bond, parameter in self.parameters:
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

                        if not self.yes_type:
                            atom.bond = 1
                        parameters_a[i] = - atom_parameter[atom.element_symbol, "A", atom.bond]

                        for j in range(i, count):
                                if i == j:
                                    distance[i, j] = atom_parameter[atom.element_symbol, "B", atom.bond]
                                else:
                                    distance[i, j] = distance[j, i] = get_distance(self.kappa,
                                                                                   molecule.atoms[i].coordinate,
                                                                                   molecule.atoms[j].coordinate)
                except KeyError:
                    print("Missing parameters for {}. element ({}{}) in {}. Program did not count with this "
                          "element.".format(atom.number, atom.element_symbol, atom.bond, name))
                    output.append((name, count, atom.element_symbol, atom.bond))
                    continue

                distance[count, :] = 1
                distance[:, count] = -1
                distance[count, count] = 0
                charges = np.linalg.solve(distance, parameters_a)
                output.append((name, count, data_from_atoms, charges))
        except KeyError:
            print("Something wrong with calculate")
        output_to_file(self.filename, self.file_parameters, output)

    def get_statistic_from_set(self):
        elements, count_element, count_all_atoms, elements_in_molecules = [], Counter(), 0, []
        for molecule in self.molecules:
            for atom in molecule.atoms:
                elements.append(atom.element_symbol)
            count_element += molecule.elements_count
            count_all_atoms += molecule.count_atoms
            for key in molecule.elements_count:
                elements_in_molecules.append(key)

            """For print statistic for singles molecule
            print("Molecule:{}{}".format(molecule.name, "".join("\n{} = {}".format(atom, count) for atom,
                                                                count in molecule.elements_count.items())))"""
        count_element_in_molecule = Counter(elements_in_molecules)
        elements_numbers = Counter(elements)
        print("Number of elements in whole set {}:".format(self.filename))
        for key in elements_numbers:
            print("{} = {}".format(key, elements_numbers[key]))
        print("Number of elements in whole set {} by bond:".format(self.filename))
        print("Element    Count in set  %in set  Found in molecules")
        for key in sorted(count_element):
            print("{0} = {1:>8} {2:12.3%} {3:>10}".format(key, count_element[key], (count_element[key]/count_all_atoms),
                                                          count_element_in_molecule[key]))

    def get_statistic_from_parameters(self):
        print("Parameters from {}:".format(self.file_parameters))
        print("Kappa = {}".format(self.kappa))
        print("Element:   A       B")
        for element, type_bond, parameter in sorted(self.parameters):
            a_parameter, b_parameter = parameter[0]
            print("{:>3}{:}: {:8} {:7}".format(element, type_bond, a_parameter, b_parameter))

    def __str__(self):
        return str("{}".format(self.molecules))


def output_to_file(filename, file_parameters, data):
    new_file = "data_from_" + filename + "_with_parameters_" + file_parameters
    with open(new_file, "w") as f:
        for name, count, atoms, charges in data:
            print("{}\n{}".format(name, count), file=f)
            for i, atom in enumerate(atoms):
                try:
                    element, number, bond = atom
                    print("{0:6d}  {1}{2:2} {3: f}".format(number, element, bond, charges[i]), file=f)
                except ValueError:
                    print("Missing parameter for this element {}({}). "
                          "Program can not calculate charges.".format(atoms, charges), file=f)
                    continue
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
    parser.add_argument("Second_Parameters_file", help="Give a file with another parameters in xml format", type=str,
                        nargs="?", default=0)
    args = parser.parse_args()
    mset = MoleculesSet()
    mset.load_parameters(args.Parameters_file)
    mset.load_from_sdf(args.SDF_file)
    mset.calculate_atom_charge()
    mset.get_statistic_from_set()
    mset.get_statistic_from_parameters()
    if args.Second_Parameters_file:
        mset.load_parameters(args.Second_Parameters_file)
        mset.calculate_atom_charge()
        mset.get_statistic_from_parameters()


if __name__ == "__main__":
    main()
