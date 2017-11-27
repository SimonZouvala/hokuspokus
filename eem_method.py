import sys
import argparse
import math
import numpy as np
import matplotlib.pyplot as plt

from collections import Counter


class Graph:
    def __init__(self, first_data):
        self.first_data = first_data

    def prepare_data(self):
        data, data_ready = [], []
        occurrence = Counter()
        for name, count, atoms, charges, parameters in self.first_data:
            for i, atom in enumerate(atoms):
                element, number, bond = atom
                element_name = element + "" + str(bond)
                occurrence[element_name] += 1
                data.append((element_name, charges[i], parameters[element, "B", bond]))
        data = sorted(data)
        start = 0
        for element in sorted(occurrence):
            prepare_data = []
            for i in range(start, occurrence[element] + start):
                element_name, charge, parameter = data[i]
                prepare_data.append((charge, parameter))
            data_ready.append((element, prepare_data))
            start += occurrence[element]
        self.data_for_draw = data_ready

    def draw_data(self):
        markers = [".", ",", "o", "v", "^", "<", ">", "p", "P", "*", "h", "+", "x", "X", "D", "d", "|", "1", "2", "3"]
        mark = 0
        for element, coordinate in sorted(self.data_for_draw):
            x, y,  = [], []
            for x_data, y_data in coordinate:
                x.append(x_data)
                y.append(y_data)
            plt.scatter(x, y, marker=markers[mark], label=element)
            mark += 1
        plt.legend(loc=1)
        plt.show()

    def draw_graph(self):
        for name, count, atoms, charges, parameters in self.first_data:
            x = np.linspace(0, count + 1, count + 1)
            for i in range(count):
                print(charges[i])
                y = x + charges[i]
        fig, ax = plt.subplots()
        ax.scatter(x, y, 'ro')
        plt.show()


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
                        elements.append(line[31:33].strip())
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

    def calculate_atom_charge(self):
        atom_parameter, output = {}, []
        self.error_molecules = 0
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
                    self.error_molecules += 1
                    print("Missing parameters for {}. element {}({}) in {}. Program did not count with this "
                          "element.".format(atom.number, atom.element_symbol, atom.bond, name))
                    output.append((name, "error", atom.element_symbol, 0, atom_parameter))
                    continue
                distance[count, :] = 1
                distance[:, count] = -1
                distance[count, count] = 0
                charges = np.linalg.solve(distance, parameters_a)
                output.append((name, count, data_from_atoms, charges, atom_parameter))
        except KeyError:
            print("Something wrong with calculate")
        self.data_for_graph = output_to_file(self.filename, self.file_parameters, output)

    def get_statistic_from_set(self):
        elements, count_element, count_all_atoms, elements_in_molecules = [], Counter(), 0, []
        for count, molecule in enumerate(self.molecules, start=1):
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
        print("Number of elements in whole set {}: {} molecules.".format(self.filename, count))
        print("Program calculated {} molecules.".format(count - self.error_molecules))
        for key in elements_numbers:
            print("{} = {}".format(key, elements_numbers[key]))
        print("Number of elements in whole set {} by bond:".format(self.filename))
        print("Element    Count in set  %in set  Found in molecules")
        for key in sorted(count_element):
            element, bond = key
            print("{0:>3}({1}) = {2:>10} {3:12.3%} {4:>10}".format(element, bond, count_element[key],
                                                                   (count_element[key]/count_all_atoms),
                                                                   count_element_in_molecule[key]))

    def get_statistic_from_parameters(self):
        print("Parameters from {}:".format(self.file_parameters))
        print("Kappa = {}".format(self.kappa))
        print("Element:   A       B")
        for element, type_bond, parameter in sorted(self.parameters):
            a_parameter, b_parameter = parameter[0]
            print("{:>3}{}: {:8} {:7}".format(element, type_bond, a_parameter, b_parameter))

    def get_data_for_graph(self):
        return self.data_for_graph

    def __str__(self):
        return str("{}".format(self.molecules))


def output_to_file(filename, file_parameters, data):
    data_for_graph = []
    new_file = "data_from_" + filename + "_with_parameters_" + file_parameters
    with open(new_file, "w") as f:
        for name, count, atoms, charges, parameters in data:
            try:
                print("{}\n{}".format(name, int(count)), file=f)
                for i, atom in enumerate(atoms):
                    element, number, bond = atom
                    print("{0:6d}  {1:>2}{2} {3: f}".format(number, element, bond, charges[i]), file=f)
                data_for_graph.append((name, count, atoms, charges, parameters))
            except ValueError:
                continue
            print("$$$$\n", file=f)
    print("Now you can find charge for each element in file {}".format(new_file))
    return data_for_graph


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
    data1 = mset.get_data_for_graph()
    mset.get_statistic_from_parameters()
    if args.Second_Parameters_file:
        mset.load_parameters(args.Second_Parameters_file)
        mset.calculate_atom_charge()
        mset.get_statistic_from_parameters()

    draw = Graph(data1)
    draw.prepare_data()
    draw.draw_data()


if __name__ == "__main__":
    main()
