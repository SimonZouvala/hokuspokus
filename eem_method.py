import sys
import argparse
import math
#import numpy


from collections import Counter, namedtuple


def get_sort_atoms_by_number_and_bond(molecules):
    number_data, atoms_sort, atoms_data, atoms_coordinate = {}, [], [], []
    for name, atoms, bonds, coordinate in molecules:
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
        number_data, atoms_data = {}, []
    #print(atoms_sort[0])
    return atoms_sort, atoms_coordinate


def load_file(file_molecules, file_parametrs):
    atoms, bonds, molecule, name, coordinate = [], [], [], "", []
    Bond = namedtuple('Bond', ['first_atom', 'second_atom', 'bond'])
    with open(file_molecules, "r") as fm:
        position, lower_position, number_name = 0, 4, 1
        data = {}
        for line in fm:
            position += 1
            if position > 3:
                data = line
                if position == 4:
                    position_bonds, position_atoms = int(data[3:6]), int(
                                                                    data[0:3])
                    atoms_rows = (position_atoms + position)
                    end_position = (atoms_rows + position_bonds + 6)
                    bonds_rows = (atoms_rows + position_bonds)
                if position <= atoms_rows and position > lower_position:
                    atoms.append(data[31:32])
                    coordinate.append((float(data[3:10]), float(data[13:20]),
                                       float(data[23:30])))

                if position > atoms_rows and position <= bonds_rows:
                    bonds.append(Bond(int(data[1:3]), int(data[3:6]),
                                 int(data[8:9])))
                if position == (end_position):
                    position_bonds, position_atoms = int(data[3:6]), int(
                                                                    data[0:3])
                    atoms_rows = (position_atoms + position)
                    bonds_rows = (atoms_rows + position_bonds)
                if "$$$$" in line:
                    molecule.append((name, atoms, bonds, coordinate))
                    number_name, end_position = position + 1, position + 4
                    lower_position, atoms, bonds = end_position, [], []
                    coordinate = []
            if position == number_name:
                name = line[0:11]
    with open(file_parametrs, "r") as fp:
        data = {}
        atoms_parametrs, parametrs = [], []
        for line in fp:
            if " Kappa=" in line:
                Kappa = float(line[line.index("Kappa") + len("Kappa=\""):-3])
            if "<Element" in line:
                atom = line[line.index("Name") + len("Name=\"")]
            if "<Bond" in line:
                parametr_A = int(line.index("A=") + len("A=\""))
                parametr_B = int(line.index("B=") + len("B=\""))
                parametrs.append((float(line[parametr_A:parametr_A + 6]),
                                float(line[parametr_B:parametr_B + 6])))
            if "</Element>" in line:
                atoms_parametrs.append((atom, parametrs))
                parametrs = []
    #print (Kappa, atoms_parametrs)
    return molecule, Kappa, atoms_parametrs


def calculate_atom_charge(atoms_sort, Kappa, atoms_parametrs, atoms_coordinate):
    atom_parametr, dx, dy, dz = {}, [], [], []
    distanse = {}
    distanse_pokus = []
    distanse_pokus2 = []
    #print (atoms_coordinate[0:2])
    for atom, parametrs in atoms_parametrs:
        #print(parametrs)
        parametr_A, parametr_B = parametrs[0]
        atom_parametr[atom, "A"] = parametr_A
        atom_parametr[atom, "B"] = parametr_B
        #print(atom_parametr[atom, "A"], atom_parametr[atom, "B"])
    for name, count, atoms in atoms_sort:
        for data in atoms:
            number, atom, bond, coordinate = data
            cx, cy, cz = coordinate
            dx.append(cx)
            dy.append(cy)
            dz.append(cz)
            #print (count, number, dx, dy, dz)
            #charge = (atom_parametr[atom, "A"] / ((-1) * atom_parametr[atom, "B"]))
            #print(charge)
        for i in range(count):
            #istanse[name, i] = atom_parametr[atoms[i][1], "B"]

            for j in range(count):
                """distanse[name, i] += Kappa * math.sqrt((
                    (dx[j] - dx[i]) ** 2) + ((dy[j] - dy[i]) ** 2)
                     + ((dz[j] - dz[i]) ** 2))
                """
                distanse[ i, j] = Kappa * math.sqrt((
                    (dx[j] - dx[i]) ** 2) + ((dy[j] - dy[i]) ** 2)
                     + ((dz[j] - dz[i]) ** 2))
                if j == count-1:
                    distanse[i, count] = -1
                if i == count-1:
                    distanse[count, j] = 1
                if i == j:
                    distanse[ i, j] = atom_parametr[atoms[i][1], "B"]
        distanse[count,count] = 0



            #print(atom_parametr[atoms[i][1],"A"],atom_parametr[atoms[i][1],"B"], distanse[name, i])
            #print(atoms[i], i,atom_parametr[atoms[i][1],"A"] / distanse[name,i] )
            #distanse_pokus2.append((i, distanse_pokus))


        print(name, distanse)
        distanse = {}
        dx, dy, dz = [], [], []

        #distanse_pokus2 = []
        #for data in atoms:
            #print(atom)

            #print (dx, dy, dz)
            #print (number, atom, bond, coordinate)
        """
        for i in range(0, count + 1):
            for j in range(0, count + 1):
                distanse_pokus.append((distanse[name,i,j]))
        """


def main():
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("Parameters", help="Give a file with molecules ",
                        type=str)
    args = parser.parse_args()
    """
    print(sys.argv[1], sys.argv[2])
    file_molecules, file_parametrs = sys.argv[1], sys.argv[2]
    molecules, Kappa, atoms_parametrs = load_file(file_molecules,
                                                  file_parametrs)
    atoms_sort, atoms_coordinate = get_sort_atoms_by_number_and_bond(molecules)
    calculate_atom_charge(atoms_sort, Kappa, atoms_parametrs, atoms_coordinate)


if __name__ == "__main__":
    main()