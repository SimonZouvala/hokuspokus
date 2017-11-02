import sys
import argparse

from collections import Counter, namedtuple


def calculate_number_of_atoms(molecules):
    counts, counts_by_bonds, count_by_one = [], Counter(), []
    for name, atoms, bonds in molecules:
        counts += [atoms for atoms in atoms]
        number_data, count_in_molecule = {}, Counter()
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
    counts = Counter(counts)
    return counts, counts_by_bonds, count_by_one


def load_file(filename):
    molecule, name = [], ""
    Bond = namedtuple('Bond', ['first_atom', 'second_atom', 'bond'])
    with open(filename, "r") as fh:
        while True:
            atoms, bonds, line = [], [], fh.readline()
            if "" == line[0:1]:
                return molecule
            name = (line[:].strip())
            for i in range(2):
                fh.readline()
            line = fh.readline()
            count_atoms, count_bonds = int(line[0:3]), int(line[3:6])
            for i in range(count_atoms):
                line = fh.readline()
                atoms.append(line[31:32])
            for i in range(count_bonds):
                line = fh.readline()
                bonds.append(Bond(int(line[1:3]), int(line[3:6]),
                                    int(line[8:9])))
            while True:
                line = fh.readline()
                if "$$$$" in line:
                    molecule.append((name, atoms, bonds))
                    break


def get_output(counts, counts_by_bonds, count_by_one):
    for name, atoms in count_by_one:
        print("Molecule: {}{}".format(name, "".join("\n{} = {}".format(atom,
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
