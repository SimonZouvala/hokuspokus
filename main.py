import argparse  # pro spouštění částí programu
import eem  # eem.py
import graph  # graph.py
import statistic  # statistic.py
import classes  # classes.py
import mgcm  # mgcm.py
import ogcm  # ogcm.py


def main():
    """ Definování všech parametrů pro spouštění různých částí programu"""
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command")
    parser_graph = subparsers.add_parser('graph')
    parser_graph.add_argument("draw_graph", type=str, nargs=2, default=0,
                              help="Give two charges results file, if you have not it, you can used CALCULATION")
    parser_calculate = subparsers.add_parser('calculation')
    parser_calculate.add_argument('calculate', type=str, nargs="?", default=0,
                                  help="Give one file with molecules (.sdf) for mgcm and "
                                       "one file with parameters for EEM (.xml)")
    parser_calculate.add_argument('--eem', action="store_true",
                                  help="Give this argument, if you want calculate with EEM")
    parser_calculate.add_argument('--mgc', action="store_true",
                                  help="Give this argument, if you want calculate with MGC")
    parser_calculate.add_argument('--ogc', action="store_true",
                                  help="Give this argument, if you want calculate with OGC")
    parser_calculate.add_argument('--parameters', type=str, help="Give this argument for parameters,"
                                                                 " if you want calculate with EEM")
    parser_calculate.add_argument('--output', type=str, help="Give a name file, for output calculate")
    parser_structure = subparsers.add_parser('structure')
    parser_structure.add_argument('--parameters', type=str, help="Give a file with parameters (EEM) (.xml)")
    parser_structure.add_argument('--molecules', type=str, help="Give a file with molecules (.sdf)")
    parser_structure.add_argument('--nobond', action="store_false",
                                  help="Give this argument, if you want not type bond.")
    args = parser.parse_args()
    try:
        if args.calculate:  # počítaní
            mset = classes.MoleculesSet()
            """eem.py"""
            if args.eem:
                set_file, para_file = args.calculate, args.parameters
                mset.load_parameters(para_file)  # načtení parametrů
                mset.load_from_sdf(set_file, args.eem, args.mgc, args.ogc)  # načtení molekuly
                cal = eem.Calculate(mset.molecules, mset.parameters)  # výpočet pomocí EEM
                """
                Parametry --eem, --mgc,  --ogc jsou boolean parametry
                """
            else:
                set_file = args.calculate
                mset.load_from_sdf(set_file, args.eem, args.mgc, args.ogc)  # načtení molekuly
            """mgcm.py"""
            if args.mgc:
                cal = mgcm.Calculate(mset.molecules, mset.periodic_table)  # výpočet pomocí MGC
            """ogcm.py"""
            if args.ogc:
                cal = ogcm.Calculate(mset.molecules)  # výpočet pomocí OGC
            if args.output:  # pro zadání uložení do souboru --output
                charge = cal
                charge.save_charges(args.output)
    except AttributeError:
        pass
    try:
        """statistic.py """
        if args.molecules:  # vypsání struktury molekul
            stat = statistic.Statistic()
            mset = classes.MoleculesSet()
            mset.load_from_sdf(args.molecules, True, False, False, args.nobond)  # True aby se chovalo jako při --eem
            stat.get_statistic_from_set(args.molecules, mset.molecules, args.nobond)  # vypsání struktury
        if args.parameters:  # vypsání struktury parametrů
            mset = classes.MoleculesSet()
            mset.load_parameters(args.parameters)  # načtení parametrů
            stat = statistic.Statistic()
            stat.get_statistic_from_parameters(args.parameters, mset.parameters)  # vypsání struktury
    except AttributeError:
        pass
    try:
        """graph.py"""
        if args.draw_graph:  # vytvoření grafu a statistických výpočtů
            graph.Graph(args.draw_graph)
    except AttributeError:
        pass


if __name__ == "__main__":
    main()
