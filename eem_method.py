import sys
import argparse
import math
import numpy as np

from collections import Counter, namedtuple


class Atom:
    def __init__(self, number, element_symbol, bond, coordinate):
        self.element_symbol = element_symbol
        self.coordinate = coordinate
        self.bond = bond
        self.number = number

    def set_bond(self, bond):
        self.bond = bond

    def __str__(self):
        return str("{} {} {} {}".format(self.number, self.element_symbol, self.bond, self.coordinate))


class Molecule:

    def __init__(self, name, atoms):
        self.name = name
        self.atoms = atoms

    def __str__(self):
        return str("Molecule name: {}\n{}".format(self.name, self.atoms))


class MoleculesSet:
    def load_from_sdf(self, filename):
        molecules = []
        with open(filename, "r") as fh:
            while True:
                number_data, atoms, elements, coordinate, line = {}, [], [], [], fh.readline()
                if "" == line[0:1]:
                    self.molecules = molecules
                    break
                name = (line[:].strip())
                for i in range(2):
                    fh.readline()
                line = fh.readline()
                count_atoms, count_bonds = int(line[0:3]), int(line[3:6])
                for i in range(1, count_atoms + 1):
                    line = fh.readline()
                    number_data[i] = 0
                    elements.append(line[31:32])
                    coordinate.append((float(line[3:10]), float(line[13:20]), float(line[23:30])))
                for i in range(count_bonds):
                    line = fh.readline()
                    first_atom, second_atom, bond = int(line[1:3]), int(line[3:6]), int(line[8:9])
                    number_data[first_atom] = max(number_data[first_atom], bond)
                    number_data[second_atom] = max(number_data[second_atom], bond)
                for key in number_data:
                    atoms.append((Atom(key, elements[key - 1], number_data[key], coordinate[key - 1])))
                while True:
                    line = fh.readline()
                    if "$$$$" in line:
                        molecules.append((Molecule(name, atoms)))
                        break

    def load_parametrs(self, filename):
        with open(filename, "r") as fp:
            parametrs = []
            for line in fp:
                if " Kappa=" in line:
                    self.kappa = float(line[line.index("Kappa") + len("Kappa=\""):-3])
                if "<Element" in line:
                    element = line[line.index("Name") + len("Name=\"")]
                    element_parametrs = []
                if "<Bond" in line:
                    parametr_a = int(line.index("A=") + len("A=\""))
                    parametr_b = int(line.index("B=") + len("B=\""))
                    element_parametrs.append((float(line[parametr_a:parametr_a + 6]),
                                              float(line[parametr_b:parametr_b + 6])))
                if "</Element>" in line:
                    parametrs.append((element, element_parametrs))
            self.parametrs = parametrs

    def __str__(self):
        return str("{}".format(self.molecules))

    def calculate_atom_charge(self):
        atom_parametr = {}
        for atom, parametr in self.parametrs:
            for bond in range(1, len(parametr) + 1):
                parametr_a, parametr_b = parametr[bond - 1]
                atom_parametr[atom, "A", bond] = parametr_a
                atom_parametr[atom, "B", bond] = parametr_b
        for number_molecule in range(len(self.molecules)):
            elements, bonds, dx, dy, dz = [], [], [], [], []
            name = self.molecules[number_molecule].name
            for number_atom in range(len(self.molecules[number_molecule].atoms)):
                count = self.molecules[number_molecule].atoms[number_atom].number
                element_symbol = self.molecules[number_molecule].atoms[number_atom].element_symbol
                bond = self.molecules[number_molecule].atoms[number_atom].bond
                coordinate = self.molecules[number_molecule].atoms[number_atom].coordinate
                cx, cy, cz = coordinate
                elements.append(element_symbol)
                bonds.append(bond)
                dx.append(cx)
                dy.append(cy)
                dz.append(cz)
            distance = np.zeros((count + 1, count+ 1))
            parametrs_a = np.zeros((count + 1, 1))
            for i in range(count):
                parametrs_a[i, 0] = atom_parametr[elements[i], "A", bonds[i]]
                for j in range(count):
                    distance[i, j] = get_distance(self.kappa, dx[j], dx[i], dy[j], dy[i], dz[j], dz[i])
                    if i == j:
                        distance[i, j] = atom_parametr[elements[i], "B", bonds[i]]
            distance[count, :] = 1
            distance[:, count] = -1
            distance[count, count] = 0
            print(name, "\n", distance, "\n", parametrs_a[:])


def get_distance(kappa, x1, x2, y1, y2, z1, z2):
    return kappa * math.sqrt(((x1 - x2) ** 2) + ((y1 - y2) ** 2) + ((z1 - z2) ** 2))


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("SDF file", help="Give a file with molecules (.sdf) ", type=str)
    parser.add_argument("Parametrs file", help="Give a file with parametrs ", type=str)
    args = parser.parse_args()
    mset = MoleculesSet()
    mset.load_parametrs(sys.argv[2])
    mset.load_from_sdf(sys.argv[1])
    mset.calculate_atom_charge()


if __name__ == "__main__":
    main()
