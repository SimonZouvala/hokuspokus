import sys
import argparse
import calculate
import graph
import statistic
import classes


from collections import Counter


def load_from_sdf(filename):
    molecules = []
    try:
        with open(filename, "r") as fh:
            while True:
                bond_data, atoms, elements, coordinate, line, elements_count = {}, [], [], [], fh.readline(), Counter()
                if "" == line[0:1]:
                    print("Load molecules from {}".format(filename))
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
                    atoms.append((classes.Atom(number, elements[number - 1], bond_data[number], coordinate[number - 1]))
                                 )
                while True:
                    line = fh.readline()
                    if "$$$$" in line:
                        molecules.append((classes.Molecule(name, count_atoms, atoms, elements_count)))
                        break
    except IOError:
        print("Wrong file for molecules set! Try another file than {}".format(filename))
        sys.exit(1)


def load_parameters(file_parameters):
    parameters_data, yes_type = [], False
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
        print("Load parameters for elements from {}".format(file_parameters))
        return parameter
    except IOError:
        print("Wrong file for parameters! Try another file than").format(file_parameters)
        sys.exit(2)


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    parser_graph = subparsers.add_parser('graph')
    parser_graph.add_argument("draw_graph", type=str, nargs=2, default=0,
                              help="Give two charges results file, if you have not it, you can used CALCULATION")
    parser_calculate = subparsers.add_parser('calculation')
    parser_calculate.add_argument('calculate', nargs=2, type=str, default=0,
                                  help="Give one file with molecules (.sdf) and one file with parameters (.xml)")
    parser_calculate.add_argument('--output', type=str, help="Give a name file, for output calculate")
    parser_statictic = subparsers.add_parser('statistic')
    parser_statictic.add_argument('--parameters', type=str, help="Give a file with parameters (.xml)")
    parser_statictic.add_argument('--molecules', type=str, help="Give a file with molecules (.sdf)")
    args = parser.parse_args()
    try:
        if args.calculate:
            set_file, para_file = args.calculate
            mset = classes.MoleculesSet(load_from_sdf(set_file))
            para = classes.Parameters(load_parameters(para_file))
            cal = calculate.Calculate(mset.molecules, para.parameters)
            if args.output:
                charge = cal
                charge.save_charges(args.output)
    except AttributeError:
        pass
    try:
        if args.molecules:
            mset = classes.MoleculesSet(load_from_sdf(args.molecules))
            stat = statistic.Statistic()
            stat.get_statistic_from_set(args.molecules, mset.molecules)
    except AttributeError:
        pass
    try:
        if args.parameters:
            para = classes.Parameters(load_parameters(args.parameters))
            stat = statistic.Statistic()
            stat.get_statistic_from_parameters(args.parameters, para.parameters)
    except AttributeError:
        pass
    try:
        if args.draw_graph:
            graph.Graph(args.draw_graph)
    except AttributeError:
        pass


if __name__ == "__main__":
    main()
