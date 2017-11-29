import sys
import main

from collections import Counter


class Load:
    def load_from_sdf(self, file_set):
        self.filename = file_set
        molecules = []
        try:
            with open(self.filename, "r") as fh:
                while True:
                    bond_data, atoms, elements, coordinate, line, elements_count = {}, [], [], [], fh.readline(), \
                                                                                   Counter()
                    if "" == line[0:1]:
                        return molecules
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
                        atoms.append((main.Atom(number, elements[number - 1], bond_data[number],
                                                coordinate[number - 1])))
                    while True:
                        line = fh.readline()
                        if "$$$$" in line:
                            molecules.append((main.Molecule(name, count_atoms, atoms, elements_count)))
                            break
        except IOError:
            print("Wrong file for molecules set! Try another file than {}".format(self.filename))
            sys.exit(1)

    def load_parameters(self, file_parameters):
        self.file_parameters, parameters_data, yes_type = file_parameters, [], False
        try:
            with open(self.file_parameters, "r") as fp:
                parameters = []
                for line in fp:
                    if " Kappa=" in line:
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
            self.parameters = parameters_data[0]
            return self.parameters
        except IOError:
            print("Wrong file for parameters! Try another file than").format(self.file_parameters)
            sys.exit(2)
