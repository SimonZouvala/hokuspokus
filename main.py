import argparse
import calculate
import graph
import output
import statistic
import loadfile


class Atom:
    def __init__(self, number, element_symbol, bond, coordinate):
        self.element_symbol = element_symbol
        self.coordinate = coordinate
        self.bond = bond
        self.number = number

    def __str__(self):
        return str("{} {} {} {}".format(self.number, self.element_symbol, self.bond, self.coordinate))


class Molecule:

    def __init__(self, name, count_atoms, atoms, elements_count):
        self.name = name
        self.count_atoms = count_atoms
        self.atoms = atoms
        self.elements_count = elements_count

    def __str__(self):
        return str("Molecule name: {}\n{}".format(self.name, self.atoms))


class MoleculesSet:
    def __init__(self, molecules):
        self.molecules = molecules

    def __str__(self):
        return str("{}".format(self.molecules))


class Parameters:
    def __init__(self, parameters):
        self.parameters = parameters


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("SDF_file", help="Give a file with molecules (.sdf) ", type=str)
    parser.add_argument("Parameters_file", help="Give a file with parameters in xml format", type=str)
    parser.add_argument("Second_Parameters_file", help="Give a file with another parameters in xml format", type=str,
                        nargs="?", default=0)
    args = parser.parse_args()
    load = loadfile.Load()
    mset = MoleculesSet(load.load_from_sdf(args.SDF_file))
    print("Load molecules from {}".format(load.filename))
    para1 = Parameters(load.load_parameters(args.Parameters_file))
    print("Load parameters for elements from {}".format(load.file_parameters))
    cal = calculate.Calculate(mset.molecules, para1.parameters)
    result, error_molecules = cal.give_result()
    out = output.Output(args.SDF_file, args.Parameters_file, result)
    file1 = out.give_name_file_with_result()
    stat = statistic.Statistic()
    stat.get_statistic_from_set(args.SDF_file, mset.molecules, error_molecules)
    stat.get_statistic_from_parameters(args.Parameters_file, para1.parameters)
    if args.Second_Parameters_file:
        para2 = Parameters(load.load_parameters(args.Second_Parameters_file))
        print("Load parameters for elements from {}".format(load.file_parameters))
        cal = calculate.Calculate(mset.molecules, para2.parameters)
        data_for_output, mset.error_molecules = cal.give_result()
        out = output.Output(args.SDF_file, args.Parameters_file, data_for_output)
        stat.get_statistic_from_parameters(args.Parameters_file, para2.parameters)
        file2 = out.give_name_file_with_result()
        draw = graph.Graph(file1, file2)
        draw.load_files()
        draw.set_graph()
        draw.graph()


if __name__ == "__main__":
    main()
