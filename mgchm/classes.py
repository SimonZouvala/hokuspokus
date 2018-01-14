import sys
import numpy as np
import re


class Atom:
    def __init__(self, number, element_symbol):
        self.element_symbol = element_symbol
        self.number = number

    def __str__(self):
        return str("{} {} ".format(self.number, self.element_symbol))


class Molecule:

    def __init__(self, name, count_atoms, atoms, count_bond_matrix, bond_matrix):
        self.name = name
        self.count_atoms = count_atoms
        self.atoms = atoms
        self.count_bond_matrix = count_bond_matrix
        self.bond_matrix = bond_matrix

    def __str__(self):
        return str("Molecule name: {}\n{}".format(self.name, self.atoms))


class MoleculesSet:
    def load_from_sdf(self, filename):
        molecules, find_elements = [], []
        try:
            with open(filename, "r") as fh:
                while True:
                    atoms, line = [], fh.readline()
                    if "" == line[0:1]:
                        self.molecules = molecules
                        self.periodic_table = get_electronegativity_from_periodic_table(set(find_elements))
                        return print("Load molecules from {}".format(filename))
                    name = (line[:].strip())
                    for i in range(2):
                        fh.readline()
                    line = fh.readline()
                    count_atoms, count_bonds = int(line[0:3]), int(line[3:6])
                    for i in range(1, count_atoms + 1):
                        line = fh.readline()
                        element = (line[31:33].strip())
                        atoms.append(Atom(i, element))
                        find_elements.append(element)
                    count_bond_matrix = np.zeros((count_atoms, count_atoms))
                    bond_matrix = np.zeros((count_atoms, count_atoms))
                    for i in range(count_bonds):
                        line = fh.readline()
                        first_atom, second_atom, bond = int(line[1:3]) - 1, int(line[3:6]) - 1, int(line[8:9])
                        bond_matrix[first_atom, second_atom] = bond_matrix[second_atom, first_atom] = bond
                        count_bond_matrix[first_atom, first_atom] += 1
                        count_bond_matrix[second_atom, second_atom] += 1
                    while True:
                        line = fh.readline()
                        if "$$$$" in line:
                            molecules.append(Molecule(name, count_atoms, atoms, count_bond_matrix, bond_matrix))
                            break
        except IOError:
            print("Wrong file for molecules set! Try another file than {}".format(filename))
            sys.exit(1)

    def __str__(self):
        return str("{}".format(self.molecules))


def get_electronegativity_from_periodic_table(elements):
    periodic_table = {}
    print(elements)
    for element in elements:
        with open("Periodic Table of Elements.csv", "r", encoding="latin-1") as ftp:
            for line in ftp:
                if ("," + element + ",") in line:
                    periodic_table[element] = float(line[[m.start() for m in re.finditer(r",", line)][10] + 1:[
                        m.start() for m in re.finditer(r",", line)][11]])
                    break
    return periodic_table
