import argparse
import calculate
import statistic
import classes


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command")
    parser_calculate = subparsers.add_parser('calculation')
    parser_calculate.add_argument('calculate', type=str, default=0,
                                  help="Give one file with molecules (.sdf) and one file with parameters (.xml)")
    parser_calculate.add_argument('--output', type=str, help="Give a name file, for output calculate")
    parser_structure = subparsers.add_parser('structure')
    parser_structure.add_argument('--molecules', type=str, help="Give a file with molecules (.sdf)")
    parser_structure.add_argument('--nobond', action="store_false", help="Give this argument, if you want"
                                                                         " not type bond.")
    args = parser.parse_args()
    print(args)

    try:
        if args.calculate:
            set_file = args.calculate
            mset = classes.MoleculesSet()
            mset.load_from_sdf(set_file)
            cal = calculate.Calculate(mset.molecules, mset.periodic_table)
            if args.output:
                charge = cal
                charge.save_charges(args.output)
    except AttributeError:
        pass
    try:
        if args.molecules:
            stat = statistic.Statistic()
            mset = classes.MoleculesSet()
            mset.load_from_sdf(args.molecules)
            stat.get_statistic_from_set(args.molecules, mset.molecules, args.nobond)
    except AttributeError:
        pass


if __name__ == "__main__":
    main()
