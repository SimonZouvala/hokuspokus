import numpy as np  # knihovna NumPy
import warnings  # pro vynechání molekuly při výskytu atomu s více než 4 vazebnými partnery
import sys  # pro ukončení práce při chybném souboru

from collections import Counter  # knihovna pro použití funkce Counter()


class Calculate:
    def __init__(self, molecules):
        self.molecules, self.output, atom_parameter, not_error_molecules = molecules, [], {}, 0
        try:
            for molecule in self.molecules:
                with warnings.catch_warnings():  # když je nějaký nedostatek při výpočtu
                    warnings.filterwarnings('error')
                    """příprava matic"""
                    tb_electronegativity, tb_hardness = molecule.tb_el, molecule.tb_hard
                    covalent_radii = molecule.tb_coval_radii
                    data_from_atoms, name, data_bond = [], molecule.name, []
                    count = molecule.count_atoms
                    degree_matrix = molecule.count_bond_matrix
                    connectivity_matrix = molecule.bond_matrix
                    identity_matrix = np.eye(degree_matrix.shape[0])
                    charges_electrons = np.zeros((degree_matrix.shape[0], 1))
                    try:
                        """výpočet S = D - A + I"""
                        simplified_matrix = degree_matrix - connectivity_matrix + identity_matrix
                        nk_electronegativity = np.linalg.solve(simplified_matrix, tb_electronegativity)
                        deviation_away_pt = nk_electronegativity-tb_electronegativity
                        denominator = 0.0
                        numerator = 0.0
                        """výpočet normalizačního faktoru Dm"""
                        for i in range(deviation_away_pt.shape[0]):
                            denominator += deviation_away_pt[i]/tb_hardness[i]
                            numerator += (deviation_away_pt[i]**2)/((tb_hardness[i]**3) * covalent_radii[i])
                        dm = denominator/numerator
                    except Warning:
                        print("Can not calculate for ", name)
                        continue
                    except Exception:
                        print("Can not calculate for ", name)
                        continue
                shift = 0
                charge_elements = Counter()
                try:
                    """Výpočet náboje jednotlivých orbitalů a výpočet náboje jednotlivých atomů"""
                    for element, number in molecule.elements_count:
                        for index in range(molecule.elements_count[element, number]):
                            """výpočet náboje orbitalu"""
                            charges_electrons[index + shift][0] = (deviation_away_pt[index + shift][0]/tb_hardness[
                                index + shift][0])-(((deviation_away_pt[index + shift][0]**2)*dm
                                                     )/((tb_hardness[index + shift][0]**3)*covalent_radii[index + shift]
                                                        ))
                            """výpočet náboje atomů"""
                            charge_elements[element, number] += float((deviation_away_pt[index + shift]/tb_hardness[
                                index + shift])-(((deviation_away_pt[index + shift]**2)*dm
                                                  )/((tb_hardness[index + shift]**3)*covalent_radii[index + shift])))
                        shift += molecule.elements_count[element, number]
                except IndexError:
                    print("Can not finished calculation for ", name)
                    continue
                not_error_molecules += 1  # součet spočítaných molekul
                for atom in molecule.atoms:  # načtení symbolu prvku a maximální vazby pro případné uložení do souboru
                    data_from_atoms.append((atom.element_symbol, atom.bond))
                self.output.append((name, count, data_from_atoms, charge_elements))
        except KeyError or IndexError:
            print("Something wrong with calculate")
            sys.exit()
        print("Program calculated {} molecules.".format(not_error_molecules))

    def save_charges(self, file):  # uložení do souboru
        data = self.output
        new_file = "result/" + file
        with open(new_file, "w") as f:
            for name, count, atoms, charges in data:
                print("{}\n{}".format(name, int(count)), file=f)
                for index, (atom, bond) in enumerate(atoms, 1):
                    print("{0:6d}  {1:>2}{2} {3: f}".format(index, atom, bond, float(charges[atom, index])), file=f)
        print("Now you can find charge for each element in file {}".format(new_file))

    def give_result(self):
        return self.output
