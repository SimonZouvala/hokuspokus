import sys
import numpy as np
import re

from collections import Counter


class Atom:
    def __init__(self, number, element_symbol, bond, coordinate=0.0):
        self.element_symbol = element_symbol
        self.coordinate = coordinate
        self.bond = bond
        self.number = number

    def __str__(self):
        return str("{} {} {} {}".format(self.number, self.element_symbol, self.bond, self.coordinate))


class Molecule:

    def __init__(self, name, count_atoms, atoms, elements_count, count_bond_matrix, bond_matrix,
                 tb_el=np.ndarray(shape=(0, 0)), tb_hard=np.ndarray(shape=(0, 0)), tb_coval_radii=np.ndarray(shape=(
                    0, 0))):
        self.name = name
        self.count_atoms = count_atoms
        self.atoms = atoms
        self.elements_count = elements_count
        self.count_bond_matrix = count_bond_matrix
        self.bond_matrix = bond_matrix
        self.tb_el = tb_el
        self.tb_hard = tb_hard
        self.tb_coval_radii = tb_coval_radii

    def __str__(self):
        return str("Molecule name: {}\n{}".format(self.name, self.atoms))


class MoleculesSet:
    def load_from_sdf(self, filename, eem, mgchm, ogchm, yes_type=True):
        molecules, find_elements = [], []
        try:
            with open(filename, "r") as fh:
                while True:
                    max_bond, atoms, elements, coordinate, line, elements_count = Counter(), [], [], [], fh.readline(),\
                                                                                  Counter()
                    shift, delete_rows, valence_state, covalent_radius = Counter(), [], {}, {}
                    atom_info = {}
                    if "" == line[0:1]:
                        self.molecules = molecules
                        if mgchm:
                            self.periodic_table = get_electronegativity_from_periodic_table(set(find_elements))
                        return print("Load molecules from {}".format(filename))
                    name = (line[:].strip())
                    for i in range(2):
                        fh.readline()
                    line = fh.readline()
                    count_atoms, count_bonds = int(line[0:3]), int(line[3:6])
                    for i in range(1, count_atoms + 1):
                        line = fh.readline()
                        elements.append(line[31:33].strip())
                        if eem:
                            coordinate.append((float(line[2:10]), float(line[12:20]), float(line[22:30])))
                        if mgchm or ogchm:
                            find_elements.append(line[31:33].strip())
                        if ogchm:
                            atom_info[i] = []
                    if ogchm:
                        size_matrix, orbital_electrons = get_orbital_electrons(elements)
                        count_bond_matrix = np.zeros((size_matrix, size_matrix))
                        bond_matrix = np.zeros((size_matrix, size_matrix))
                        table_electronegativity = np.zeros((size_matrix, 1))
                        table_hardness = np.zeros((size_matrix, 1))
                        matrix_covalent_radii = np.zeros((size_matrix, 1))
                    if mgchm:
                        count_bond_matrix = np.zeros((count_atoms, count_atoms))
                        bond_matrix = np.zeros((count_atoms, count_atoms))
                    for i in range(count_bonds):
                        line = fh.readline()
                        first_atom, second_atom, bond = int(line[1:3]), int(line[3:6]), int(line[8:9])
                        max_bond[first_atom] = max(max_bond[first_atom], bond)
                        max_bond[second_atom] = max(max_bond[second_atom], bond)
                        if mgchm:
                            index1, index2 = first_atom - 1, second_atom - 1
                            bond_matrix[index1, index2] = bond_matrix[index2, index1] = bond
                            count_bond_matrix[index1, index1] += bond
                            count_bond_matrix[index2, index2] += bond
                        if ogchm:
                            info_first, info_second = list(atom_info[first_atom]), list(atom_info[second_atom])
                            info_first.append(bond)
                            info_second.append(bond)
                            atom_info[first_atom] = info_first
                            atom_info[second_atom] = info_second
                            orbital_electrons[0] = 0
                            position1 = position2 = 0
                            for valence in range(0, first_atom):
                                position1 += orbital_electrons[valence]
                            for valence in range(0, second_atom):
                                position2 += orbital_electrons[valence]
                        try:
                            for x in range(bond):
                                bond_matrix[position1 + shift[first_atom]][position2 + shift[second_atom]] = \
                                    bond_matrix[position2 + shift[second_atom]][position1 + shift[first_atom]] = 1
                                shift[first_atom] += 1
                                shift[second_atom] += 1
                        except IndexError:
                            print("Can not prepare data from ", name)
                            break
                    try:
                        if ogchm:
                            orbital_electrons, max_bond, table_electronegativity, table_hardness,\
                            matrix_covalent_radii, valence_state, \
                            delete_rows = get_valence_state_and_prepare_table_values(atom_info, elements,
                                                                                     orbital_electrons, max_bond,
                                                                                     table_electronegativity,
                                                                                     table_hardness,
                                                                                     matrix_covalent_radii)
                            deleted = 0
                            for row in sorted(set(delete_rows)):
                                table_electronegativity = np.delete(table_electronegativity, row - deleted, 0)
                                table_hardness = np.delete(table_hardness, row - deleted, 0)
                                matrix_covalent_radii = np.delete(matrix_covalent_radii, row - deleted, 0)
                                for columns_or_row in [0, 1]:
                                    bond_matrix = np.delete(bond_matrix, row - deleted, columns_or_row)
                                    count_bond_matrix = np.delete(count_bond_matrix, row - deleted, columns_or_row)
                                deleted += 1
                            position = 0
                            for valence_electron in orbital_electrons:
                                if orbital_electrons[valence_electron] >= (shift[valence_electron] + 1):
                                    while orbital_electrons[valence_electron] > 4:
                                        orbital_electrons[valence_electron] -= 1
                                for electron1 in range(orbital_electrons[valence_electron]):
                                    for electron2 in range(orbital_electrons[valence_electron]):
                                        if electron1 != electron2:
                                            bond_matrix[electron1 + position, electron2 + position] = 1
                                position += orbital_electrons[valence_electron]
                            for index in range(count_bond_matrix.shape[0]):
                                count_bond_matrix[index, index] = int(sum(bond_matrix[index, :]))
                    except IndexError:
                        print("Can not calculate for ", name)
                        while True:
                            line = fh.readline()
                            if "$$$$" in line:
                                break
                    for number in sorted(max_bond):
                        if eem:
                            if yes_type:
                                elements_count[(elements[number - 1], max_bond[number])] += 1
                            else:
                                elements_count[elements[number - 1]] += 1
                            atoms.append(Atom(number, elements[number - 1], max_bond[number], coordinate[number - 1]))
                        else:
                            if ogchm:
                                elements_count[elements[number - 1], number] = len(valence_state[number])
                            atoms.append(Atom(number, elements[number - 1], max_bond[number]))
                    while True:
                        line = fh.readline()
                        if "$$$$" in line:
                            if eem:
                                count_bond_matrix = bond_matrix = False
                            if mgchm:
                                elements_count = False
                            if ogchm:
                                molecules.append(Molecule(name, count_atoms, atoms, elements_count, count_bond_matrix,
                                                          bond_matrix, table_electronegativity, table_hardness,
                                                          matrix_covalent_radii))
                            else:
                                molecules.append(Molecule(name, count_atoms, atoms, elements_count, count_bond_matrix,
                                                          bond_matrix))
                            break
        except IOError:
            print("Wrong file for molecules set! Try another file than {}".format(filename))
            sys.exit(1)

    def load_parameters(self, file_parameters):
        parameters_data, yes_type, element, kappa = [], False, "", 0.0
        try:
            with open(file_parameters, "r") as fp:
                parameters = []
                for line in fp:
                    if "Kappa=" in line:
                        kappa = float(line[line.index("Kappa") + len("Kappa=\""):-3])
                    if "<Element" in line:
                        element = line[line.index("Name") + len("Name=\""):-3].strip()
                    if "<Bond" in line:
                        element_parameters = []
                        if "Type" in line:
                            yes_type = True
                            type_bond = int(line[line.index("Type") + len("Type=\"")])
                        else:
                            type_bond = 1
                        parameter_a = int(line.index("A=") + len("A=\""))
                        parameter_b = int(line.index("B=") + len("B=\""))
                        element_parameters.append((float(line[parameter_a:parameter_a + 6]),
                                                   float(line[parameter_b:parameter_b + 6])))
                        parameters.append((element, type_bond, element_parameters))
            parameters_data.append((kappa, yes_type, parameters))
            parameter = parameters_data[0]
            self.parameters = parameter
            print("Load parameters for elements from {}".format(file_parameters))
        except IOError:
            print("Wrong file for parameters! Try another file than").format(file_parameters)
            sys.exit(2)

    def __str__(self):
        return str("{}".format(self.molecules))


def get_electronegativity_from_periodic_table(elements):
    periodic_table = {}
    for element in elements:
        with open("data/Periodic Table of Elements.csv", "r", encoding="latin-1") as ftp:
            for line in ftp:
                if ("," + element + ",") in line:
                    periodic_table[element] = float(line[[m.start() for m in re.finditer(r",", line)][10] + 1:[
                        m.start() for m in re.finditer(r",", line)][11]])
                    break
    return periodic_table


def get_orbital_electrons(elements):
    matrix, orbital_electrons, ionization_potential = 0, Counter(), {}
    for x, element in enumerate(elements):
        x += 1
        with open("data/Periodic Table of Elements.csv", "r", encoding="latin-1") as ftp:
            for line in ftp:
                if ("," + element + ",") in line:
                    number = int(line[0:[m.start() for m in re.finditer(r",", line)][0]])
                    if number < 3:
                        valence_number = int(line[[m.start() for m in re.finditer(r",", line)][-3] + 3:[
                            m.start() for m in re.finditer(r",", line)][-2]])
                    else:
                        valence_number = 0
                        valence_text = line[[m.start() for m in re.finditer(r",", line)][-3] + 5:[
                            m.start() for m in re.finditer(r",", line)][-2]]
                        if len(valence_text) > 8:
                            if len(valence_text) > 10:
                                valence_text = valence_text[-8:]
                                if len(valence_text) > 8:
                                    valence_text = valence_text[-1]
                                    valence_number = int(valence_text[-1])
                            else:
                                valence_text = valence_text[-1]
                                valence_number = int(valence_text[-1])
                                break
                        maximum = int((len(valence_text)) / 4)
                        for i in range(maximum):
                            number = valence_text[3:5].strip()
                            valence_number += int(number)
                            valence_text = valence_text[3 + len(number):]
                    orbital_electrons[x] = valence_number
                    matrix += valence_number
                    break
    return matrix, orbital_electrons


def get_electronnegativity_and_hardness(state_text, element):
    table_electronegativity, table_hardness = Counter(), Counter()
    with open("data/electronegativity_hardness.csv", "r", encoding="latin-1") as feh:
        for line in feh:
            if (element + ",") in line:
                if set(state_text) == set((line[[m.start() for m in re.finditer(r",", line)][0] + 1:[
                            m.start() for m in re.finditer(r",", line)][1]]).split(" ")):
                    table_electronegativity[element, line[[m.start() for m in re.finditer(r",", line)][1] + 1:[
                                m.start() for m in re.finditer(r",", line)][2]]] = float(
                        line[[m.start() for m in re.finditer(r",", line)][2] + 1:[
                                m.start() for m in re.finditer(r",", line)][3]])
                    table_hardness[element, line[[m.start() for m in re.finditer(r",", line)][1] + 1:[
                                m.start() for m in re.finditer(r",", line)][2]]] = float(
                        line[[m.start() for m in re.finditer(r",", line)][3] + 1:])
    return table_electronegativity, table_hardness


def get_covalent_radii(element, bond):
    covalent_radius = 0
    with open("data/Covalent radii.csv", "r", encoding="latin-1") as fcr:
        for line in fcr:
            if (element + "," + str(bond)) in line:
                covalent_radius = float(line[[m.start() for m in re.finditer(r",", line)][1] + 1:])
    return covalent_radius


def get_valence_state_and_prepare_table_values(bond_info, elements, orbital_electrons, max_bond,
                                               table_electronegativity, table_hardness, matrix_covalent_radii):
    covalent_radius, valence_state, delete_rows = {}, {}, []
    orbital_electrons[0] = 0
    for number_element in bond_info:
        position = 0
        for valence in range(0, number_element):
            position += orbital_electrons[valence]
        covalent_radius[elements[number_element - 1], max(bond_info[number_element])] = get_covalent_radii(
            elements[number_element - 1], max(bond_info[number_element]))
        state_text, filled_type_state, non_binding_pair, sigma = [], 1, 0, 0
        valence_electrons = orbital_electrons[number_element]
        prepare = valence_electrons - sum(bond_info[number_element])
        long_state, index_state, text, filled = valence_electrons, 4 - (min(4, valence_electrons) - 1), "", 0
        if prepare/2 >= 1:
            long_state = int(valence_electrons - (prepare/2))
        if sum(bond_info[number_element]) - len(bond_info[number_element]) > 0:
            index_state = sum(bond_info[number_element]) - len(bond_info[number_element]) + 1
        electron_states = ["s", "di", "tr", "te"]
        divisor = 2
        for bond in bond_info[number_element]:
            if bond > 1:
                for x in range(bond - 1):
                    state_text.append("pi")
                    filled_type_state += 1
                    filled += 1
                state_text.append(electron_states[-index_state])
                filled_type_state += 1
                filled += 1
                divisor = 2
            else:
                state_text.append(electron_states[-index_state])
                filled += 1
        if ((valence_electrons - filled) / divisor) >= 1:
            for y in range((valence_electrons - filled) - (long_state - filled)):
                state_text.append(electron_states[-index_state] + "2")
                non_binding_pair += 1
        electronegativity, hardness = get_electronnegativity_and_hardness(state_text, elements[number_element - 1])
        valence_state[number_element] = state_text
        dive = 1
        for move, state in enumerate(state_text):
            table_electronegativity[position + move][0] = electronegativity[
                elements[number_element - 1], state]
            table_hardness[position + move][0] = hardness[elements[number_element - 1], state]
            matrix_covalent_radii[position + move][0] = covalent_radius[elements[number_element - 1],
                                                                        max(bond_info[number_element])]
            if "2" in state:
                if move <= 3:
                    dive += 1
                delete_rows.append(int(position + len(state_text) - 2 + dive))
    return orbital_electrons, max_bond, table_electronegativity, table_hardness, matrix_covalent_radii, valence_state,\
           delete_rows
