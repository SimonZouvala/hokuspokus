from collections import Counter


class Statistic:
    def get_statistic_from_set(self, file_set, molecules, type_bond=True):
        elements, count_element, count_all_atoms, elements_in_molecules = [], Counter(), 0, []
        for count, molecule in enumerate(molecules, start=1):
            for atom in molecule.atoms:
                elements.append(atom.element_symbol)
            count_element += molecule.elements_count
            count_all_atoms += molecule.count_atoms
            for key in molecule.elements_count:
                elements_in_molecules.append(key)
        count_element_in_molecule = Counter(elements_in_molecules)
        print("Number of elements in whole set {}: {} molecules.".format(file_set, count))

        """----For print number of elements in whole set not by bond----
        elements_numbers = Counter(elements)
        for key in elements_numbers:
            print("{} = {}".format(key, elements_numbers[key]))
        """
        print("Element    Count in set  %in set  Found in molecules")
        if type_bond:
            for key in sorted(count_element):
                element, bond = key
                print("{0:>3}({1}) = {2:>10} {3:12.3%} {4:>10}".format(element, bond, count_element[key],
                                                                       (count_element[key] / count_all_atoms),
                                                                       count_element_in_molecule[key]))
        else:
            for element in sorted(count_element):
                print("{0:>3} = {1:>10} {2:12.3%} {3:>10}".format(element, count_element[element],
                                                                  (count_element[element] / count_all_atoms),
                                                                  count_element_in_molecule[element]))

    def get_statistic_from_parameters(self, file_parameters, parameters):
        kappa, yes_type, parameters = parameters
        print("Parameters from {}:".format(file_parameters))
        print("Kappa = {}".format(kappa))
        print("Element:   A       B")
        for element, type_bond, parameter in sorted(parameters):
            a_parameter, b_parameter = parameter[0]
            print("{:>3}{}: {:8} {:7}".format(element, type_bond, a_parameter, b_parameter))
