import argparse
import calculate
import graph
import statistic
import classes


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
    parser_structure = subparsers.add_parser('structure')
    parser_structure.add_argument('--parameters', type=str, help="Give a file with parameters (.xml)")
    parser_structure.add_argument('--molecules', type=str, help="Give a file with molecules (.sdf)")
    parser_structure.add_argument('--bond', action="store_false", help="Give this argument, if you want not type bond.")
    args = parser.parse_args()
    try:
        if args.calculate:
            set_file, para_file = args.calculate
            mset = classes.MoleculesSet()
            mset.load_from_sdf(set_file)
            mset.load_parameters(para_file)
            cal = calculate.Calculate(mset.molecules, mset.parameters)
            if args.output:
                charge = cal
                charge.save_charges(args.output)
    except AttributeError:
        pass
    try:
        if args.molecules:
            stat = statistic.Statistic()
            mset = classes.MoleculesSet()
            mset.load_from_sdf(args.molecules, args.bond)
            stat.get_statistic_from_set(args.molecules, mset.molecules, args.bond)
    except AttributeError:
        pass
    try:
        if args.parameters:
            mset = classes.MoleculesSet()
            mset.load_from_sdf(args.parameters)
            stat = statistic.Statistic()
            stat.get_statistic_from_parameters(args.parameters, mset.parameters)
    except AttributeError:
        pass
    try:
        if args.draw_graph:
            graph.Graph(args.draw_graph)
    except AttributeError:
        pass


if __name__ == "__main__":
    main()
