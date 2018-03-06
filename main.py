import argparse
import eem
import graph
import statistic
import classes
import mgchm
import ogchm


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command")
    parser_graph = subparsers.add_parser('graph')
    parser_graph.add_argument("draw_graph", type=str, nargs=2, default=0,
                              help="Give two charges results file, if you have not it, you can used CALCULATION")
    parser_calculate = subparsers.add_parser('calculation')
    parser_calculate.add_argument('calculate', type=str, nargs='?', default=0,
                                  help="Give one file with molecules (.sdf) for MGCHM and "
                                       "one file with parameters for EEM (.xml)")
    parser_calculate.add_argument('--eem', action="store_true",
                                  help="Give this argument, if you want calculate with EEM")
    parser_calculate.add_argument('--mgchm', action="store_true",
                                  help="Give this argument, if you want calculate with MGCHM")
    parser_calculate.add_argument('--ogchm', action="store_true",
                                  help="Give this argument, if you want calculate with OGCHM")
    parser_calculate.add_argument('--output', type=str, help="Give a name file, for output calculate")
    parser_structure = subparsers.add_parser('structure')
    parser_structure.add_argument('--parameters', type=str, help="Give a file with parameters (EEM) (.xml)")
    parser_structure.add_argument('--molecules', type=str, help="Give a file with molecules (.sdf)")
    parser_structure.add_argument('--nobond', action="store_false",
                                  help="Give this argument, if you want not type bond.")
    args = parser.parse_args()
    print(args)
    try:
        if args.calculate:
            mset = classes.MoleculesSet()
            if args.eem:
                set_file, para_file = args.calculate
                mset.load_parameters(para_file)
                mset.load_from_sdf(set_file, args.eem, args.mgchm, args.ogchm)
                cal = eem.Calculate(mset.molecules, mset.parameters)
            else:
                set_file = args.calculate
                mset.load_from_sdf(set_file, args.eem, args.mgchm, args.ogchm)
            if args.mgchm:
                cal = mgchm.Calculate(mset.molecules, mset.periodic_table)
            if args.ogchm:
                cal = ogchm.Calculate(mset.molecules, mset.tb_el, mset.tb_hard)
            if args.output:
                charge = cal
                charge.save_charges(args.output)
    except AttributeError:
        pass
    try:
        if args.molecules:
            stat = statistic.Statistic()
            mset = classes.MoleculesSet()
            mset.load_from_sdf(args.molecules, True, False, False, args.nobond)
            stat.get_statistic_from_set(args.molecules, mset.molecules, args.nobond)
        if args.parameters:
            mset = classes.MoleculesSet()
            mset.load_parameters(args.parameters)
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
