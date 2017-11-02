import sys
import argparse
import math
#import numpy


from collections import Counter, namedtuple


def get_sort_atoms_by_number_and_bond(molecules):
    atoms_sort, atoms_coordinate = [], []
    for name, atoms, bonds, coordinate in molecules:
        number_data, atoms_data = {}, []
        for i in range(1, len(atoms) + 1):
            number_data[i] = 0
        for i in range(len(bonds)):
            first_atom, second_atom, bond = bonds[i]
            number_data[first_atom] = max(number_data[first_atom], bond)
            number_data[second_atom] = max(number_data[second_atom], bond)
        for key in number_data:
            atoms_data.append((key, atoms[key - 1], number_data[key],
                               coordinate[key - 1]))
            atoms_coordinate.append((key, coordinate[key - 1]))
        atoms_sort.append((name, len(atoms), atoms_data))
    print(atoms_coordinate)
    return atoms_sort, atoms_coordinate


def load_file(file_molecules, file_parametrs):
    molecule, name, atoms_parametrs = [], "", []
    Bond = namedtuple('Bond', ['first_atom', 'second_atom', 'bond'])
    with open(file_parametrs, "r") as fp:
        for line in fp:
            if " Kappa=" in line:
                Kappa = float(line[line.index("Kappa") + len("Kappa=\""):-3])
            if "<Element" in line:
                atom = line[line.index("Name") + len("Name=\"")]
                parametrs = []
            if "<Bond" in line:
                parametr_A = int(line.index("A=") + len("A=\""))
                parametr_B = int(line.index("B=") + len("B=\""))
                parametrs.append((float(line[parametr_A:parametr_A + 6]),
                                float(line[parametr_B:parametr_B + 6])))
            if "</Element>" in line:
                atoms_parametrs.append((atom, parametrs))
    with open(file_molecules, "r") as fh:
        while True:
            atoms, bonds, coordinate, line = [], [], [], fh.readline()
            if "" == line[0:1]:
                return molecule, Kappa, atoms_parametrs
            name = (line[:].strip())
            for i in range(2):
                fh.readline()
            line = fh.readline()
            count_atoms, count_bonds = int(line[0:3]), int(line[3:6])
            for i in range(count_atoms):
                line = fh.readline()
                atoms.append(line[31:32])
                coordinate.append((float(line[3:10]), float(line[13:20]),
                                       float(line[23:30])))
            for i in range(count_bonds):
                line = fh.readline()
                bonds.append(Bond(int(line[1:3]), int(line[3:6]),
                                  int(line[8:9])))
            while True:
                line = fh.readline()
                if "$$$$" in line:
                    molecule.append((name, atoms, bonds, coordinate))
                    break


def calculate_atom_charge(atoms_sort, Kappa, atoms_parametrs, atoms_coordinate):
    atom_parametr, dx, dy, dz = {}, [], [], []
    distanse = {}
    for atom, parametr in atoms_parametrs:
        print (parametr)
        parametr_A, parametr_B = parametr[0]
        atom_parametr[atom, "A"] = parametr_A
        atom_parametr[atom, "B"] = parametr_B
    for name, count, atoms in atoms_sort:
        for data in atoms:
            number, atom, bond, coordinate = data
            cx, cy, cz = coordinate
            dx.append(cx)
            dy.append(cy)
            dz.append(cz)
        for i in range(count):
            for j in range(count):
                """distanse[name, i] += Kappa * math.sqrt((
                    (dx[j] - dx[i]) ** 2) + ((dy[j] - dy[i]) ** 2)
                     + ((dz[j] - dz[i]) ** 2))
                """
                distanse[i, j] = Kappa * math.sqrt((
                    (dx[j] - dx[i]) ** 2) + ((dy[j] - dy[i]) ** 2)
                     + ((dz[j] - dz[i]) ** 2))
                if j == count - 1:
                    distanse[i, count] = -1
                if i == count - 1:
                    distanse[count, j] = 1
                if i == j:
                    distanse[i, j] = atom_parametr[atoms[i][1], "B"]
        distanse[count, count] = 0
        print(name, distanse)
        distanse = {}
        dx, dy, dz = [], [], []


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("SDF file", help="Give a file with molecules (.sdf) ",
                        type=str)
    parser.add_argument("Parametrs file", help="Give a file with parametr ",
                        type=str)
    args = parser.parse_args()
    file_molecules, file_parametrs = sys.argv[1], sys.argv[2]
    molecules, Kappa, atoms_parametrs = load_file(file_molecules,
                                                  file_parametrs)
    atoms_sort, atoms_coordinate = get_sort_atoms_by_number_and_bond(molecules)
    calculate_atom_charge(atoms_sort, Kappa, atoms_parametrs, atoms_coordinate)


if __name__ == "__main__":
    main()