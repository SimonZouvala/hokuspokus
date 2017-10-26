import sys
import argparse

from collections import Counter, namedtuple


def calculate_number_of_atoms(molecules):
    counts, number_data, counts_by_bonds, = [], {}, Counter()
    count_in_molecule, count_by_one = Counter(), []
    for name, atoms, bonds in molecules:
        counts += [atoms for atoms in atoms]
        for i in range(1, len(atoms) + 1):
            number_data[i] = 0
        for i in range(len(bonds)):
            first_atom, second_atom, bond = bonds[i]
            number_data[first_atom] = max(number_data[first_atom], bond)
            number_data[second_atom] = max(number_data[second_atom], bond)
        for key in number_data:
            counts_by_bonds[(atoms[key - 1], number_data[key])] += 1
            count_in_molecule[(atoms[key - 1], number_data[key])] += 1
        count_by_one.append((name, count_in_molecule))
        number_data, count_in_molecule = {}, Counter()
    counts = Counter(counts)
    return counts, counts_by_bonds, count_by_one


def load_file(filename):
    atoms, bonds, molecule, name = [], [], [], ""
    Bond = namedtuple('Bond', ['first_atom', 'second_atom', 'bond'])
    with open(filename, "r") as fh:
        position, lower_position, number_name = 0, 4, 1
        data = {}
        for line in fh:
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
                if position > atoms_rows and position <= bonds_rows:
                    bonds.append(Bond(int(data[1:3]), int(data[3:6]),
                                 int(data[8:9])))
                if position == (end_position):
                    position_bonds, position_atoms = int(data[3:6]), int(
                                                                    data[0:3])
                    atoms_rows = (position_atoms + position)
                    bonds_rows = (atoms_rows + position_bonds)
                if "$$$$" in line:
                    molecule.append((name, atoms, bonds))
                    number_name, end_position = position + 1, position + 4
                    lower_position, atoms, bonds = end_position, [], []
            if position == number_name:
                name = line[0:11]
    return molecule


def get_output(counts, counts_by_bonds, count_by_one):
    for name, atoms in count_by_one:
        print("Molecule: {}{}".format(name, "".join("{} = {}\n".format(atom,
                                      count) for atom, count in atoms.items())))
    print ("Number of atoms in the whole set:")
    for key in counts:
        print("{} = {}".format(key, counts[key]))
    print ("Number of atoms by bonds:")
    for atom in sorted(counts_by_bonds):
        if counts_by_bonds[atom] > 0:
            print("{} = {}".format(atom, counts_by_bonds[atom]))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("Parameters", help="Give a file with molecules",
                        type=str)
    args = parser.parse_args()
    filename = sys.argv[1]
    molecules = load_file(filename)
    counts, counts_by_bonds, count_by_one = calculate_number_of_atoms(molecules
                                                                      )
    get_output(counts, counts_by_bonds, count_by_one)

if __name__ == "__main__":
    main()
