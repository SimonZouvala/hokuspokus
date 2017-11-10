import sys
import argparse
import math
import numpy as np

from collections import Counter, namedtuple


class Parametrs:

    def set_kappa(self, kappa):
        self.kappa = kappa

    def __init__(self, atom, parametr_a, parametr_b):
        self.atom = atom
        self.parametr_a = parametr_a
        self.parametr_b = parametr_b

    pass


class Atom:
    def __init__(self, i, element_symbol, bond, coordinate):
        self.element_symbol = element_symbol
        self.coordinate = coordinate
        self.bond = bond
        self.i = i

    def set_bond(self, bond):
        self.bond = bond

    def __str__(self):
        return str("{} {} {} {}".format(self.i, self.element_symbol, self.bond, self.coordinate))


class Molecule:
    #atoms: list[Atom]

    def __init__(self, name):
        self.name = name

    def get_sort_atoms_by_number_and_bond(self, data):
        atoms_sort, atoms_coordinate = [], []
        for name, atoms, bonds, coordinate in data:
            number_data, atoms_data = {}, []
            for i in range(1, len(atoms) + 1):


                atoms_coordinate.append((key, coordinate[key - 1]))
            atoms_sort.append((name, len(atoms), atoms_data))
        return atoms_sort, atoms_coordinate
    pass


class MoleculesSet:
    #molecules:  list[Molecule]

    def load_from_sdf(self, filename):
        Bond = namedtuple('Bond', ['first_atom', 'second_atom', 'bond'])
        i = 0
        with open(filename, "r") as fh:
            while True:
                number_data, atoms, bonds, coordinate, line = {}, [], [], [], fh.readline()
                if "" == line[0:1]:
                    break
                name = (line[:].strip())
                for i in range(2):
                    fh.readline()
                line = fh.readline()
                count_atoms, count_bonds = int(line[0:3]), int(line[3:6])
                for i in range(count_atoms):
                    line = fh.readline()
                    i += 1
                    number_data[i] = 0
                    atoms.append(line[31:32])
                    coordinate.append((float(line[3:10]), float(line[13:20]),
                                float(line[23:30])))
                for i in range(count_bonds):
                    line = fh.readline()
                    first_atom, second_atom, bond = int(line[1:3]), int(line[3:6]), int(line[8:9])
                    number_data[first_atom] = max(number_data[first_atom], bond)
                    number_data[second_atom] = max(number_data[second_atom], bond)
                for key in number_data:
                    atom = Atom(key, atoms[key - 1], number_data[key], coordinate[key - 1])
                    print(atom)
                while True:
                    line = fh.readline()
                    if "$$$$" in line:
                        molecule = Molecule(name)
                        i = 0
                        break

    def load_parametrs(self, filename):
        with open(filename, "r") as fp:
            for line in fp:
                if " Kappa=" in line:
                    kappa = float(line[line.index("Kappa") + len("Kappa=\""):-3])
                if "<Element" in line:
                    element = line[line.index("Name") + len("Name=\"")]
                if "<Bond" in line:
                    parametr_a = int(line.index("A=") + len("A=\""))
                    parametr_b = int(line.index("B=") + len("B=\""))
                    a = line[parametr_a + 1:parametr_a + 5]
                    b = line[parametr_b + 1:parametr_b + 5]
                if "</Element>" in line:
                    parametrs = Parametrs(element, a, b)
                    parametrs.set_kappa(kappa)
    pass


def calculate_atom_charge(atoms_sort, Kappa, atoms_parametrs, atoms_coordinate):
    atom_parametr, dx, dy, dz = {}, [], [], []
    for atom, parametr in atoms_parametrs:
        parametr_A, parametr_B = parametr[0]
        atom_parametr[atom, "A"] = parametr_A
        atom_parametr[atom, "B"] = parametr_B
    for name, count, atoms in atoms_sort:
        element = []
        identity_matrix = np.eye(count+1)
        distance = np.zeros((count + 1, count + 1))
        parametrsA = np.zeros((count+1, 1))
        for data in atoms:
            number, atom, bond, coordinate = data
            cx, cy, cz = coordinate
            element.append(atom)
            dx.append(cx)
            dy.append(cy)
            dz.append(cz)
        for i in range(len(element)):
            parametrsA[i, 0] = atom_parametr[element[i], "A"]
        for i in range(count):
            for j in range(count):
                distance[i, j] = Kappa * math.sqrt((
                    (dx[j] - dx[i]) ** 2) + ((dy[j] - dy[i]) ** 2)
                     + ((dz[j] - dz[i]) ** 2))
                if i == j:
                    distance[i, j] = atom_parametr[atoms[i - 1][1], "B"]
        distance[count, :] = 1
        distance[:, count] = -1
        distance[count, count] = 0
        print(name, distance)
        dx, dy, dz = [], [], []


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("SDF file", help="Give a file with molecules (.sdf) ",
                        type=str)
    parser.add_argument("Parametrs file", help="Give a file with parametr ",
                        type=str)
    args = parser.parse_args()
    mset = MoleculesSet()
    mset.load_parametrs(sys.argv[2])
    mset.load_from_sdf(sys.argv[1])


"""
    file_molecules, file_parametrs = sys.argv[1], sys.argv[2]
    molecules, Kappa, atoms_parametrs = load_file(file_molecules,
                                                  file_parametrs)
    atoms_sort, atoms_coordinate = get_sort_atoms_by_number_and_bond(molecules)
    calculate_atom_charge(atoms_sort, Kappa, atoms_parametrs, atoms_coordinate)
"""

if __name__ == "__main__":
    main()
