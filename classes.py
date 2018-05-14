import sys  # pro ukončení práce při chybném souboru
import numpy as np  # NumPy knihovna
import re  # pro pomoc vyhledávání znaku v souboru

from collections import Counter  # knihovna pro použití funkce Counter()


"""TŘÍDA Atom"""


class Atom:
    def __init__(self, number, element_symbol, bond, coordinate=0.0):
        self.element_symbol = element_symbol  # symbol prvku
        self.coordinate = coordinate  # souřadnice atomu
        self.bond = bond  # vazba
        self.number = number  # číslo atomu

    def __str__(self):
        return str("{} {} {} {}".format(self.number, self.element_symbol, self.bond, self.coordinate))


"""TŘÍDA Molecule"""


class Molecule:

    def __init__(self, name, count_atoms, atoms, elements_count, count_bond_matrix, bond_matrix,
                 tb_el=np.ndarray(shape=(0, 0)), tb_hard=np.ndarray(shape=(0, 0)), tb_coval_radii=np.ndarray(shape=(
                    0, 0))):
        self.name = name  # název molekuly
        self.count_atoms = count_atoms  # počet atomů v molekule
        self.atoms = atoms  # atomy v molekule, které odkazují na třídu Atom
        self.elements_count = elements_count  # počet atomů o maximální vazbě (EEM) nebo uložení valenční vrstvy (OGC)
        self.count_bond_matrix = count_bond_matrix  # matice stupně (MGC, OGC)
        self.bond_matrix = bond_matrix  # matice sousednosti (MGC, OGC)
        self.tb_el = tb_el  # vektor tabulkových hodnot Elektronegativity (OGC)
        self.tb_hard = tb_hard  # vektro tabulkových hodnot tvrdosti (OGC)
        self.tb_coval_radii = tb_coval_radii  # vektor tabulkových hodnot kovalentního poloměru (OGC)

    def __str__(self):
        return str("Molecule name: {}\n{}".format(self.name, self.atoms))


"""NAČTENÍ MOLEKULY DO TŘÍDY MoleculesSet"""


class MoleculesSet:
    def load_from_sdf(self, filename, eem, mgc, ogc, yes_type=True):
        """
        molecules slouží pro uložení informací o atomech v molekule
        find_element pro zjistění elektronegativit atomů, které jsou v sadě molekul při výpočetu MGC
        """
        molecules, find_elements = [], []
        try:
            with open(filename, "r") as fh:
                while True:
                    max_bond, atoms, elements, coordinate, line, elements_count = Counter(), [], [], [], fh.readline(),\
                                                                                  Counter()
                    """ 
                    max_bond pro zjištění maximální vazby u atomu
                    atoms pro uložení třídy Atom pro jednotlivé molekuly
                    elemnts pro uložení jednotlivých prvků
                    coordinate pro uložení 3D souřadnic při počítání s EEM
                    elements_count pro získání počtu atomů s jednotlivýma vazbama v molekule
                    shift pro počítání posunu na další elektron, který ještě není ve vazbě, u jednotlivých atomů
                    delete_rows pro zjištění, který řádky a sloupce s ebudou odstraňovat z matic
                    valence_state pro uložení valenčního stavu atomu
                    atom_info pro uložení všech vazeb jednotlivých atomů
                    """
                    shift, delete_rows, valence_state = Counter(), [], {}
                    atom_info = {}
                    if "" == line[0:1]:  # Na konci souboru, se uloží molekuly do třídy molecules
                        self.molecules = molecules
                        if mgc:
                            """Při mětodě MGC se zjistí elektronegativita pro všechny atomy, které jsou v molekulách"""
                            self.periodic_table = get_electronegativity_from_periodic_table(set(find_elements))
                        return print("Load molecules from {}".format(filename))
                    name = (line[:].strip())
                    for i in range(2):  # Nepotřebné řádky
                        fh.readline()
                    line = fh.readline()
                    count_atoms, count_bonds = int(line[0:3]), int(line[3:6])  # Počet atomů a vazeb v molekule
                    for i in range(1, count_atoms + 1):  # Načtení atomů, souřadnic při EEM
                        line = fh.readline()
                        elements.append(line[31:33].strip())
                        if eem:
                            coordinate.append((float(line[2:10]), float(line[12:20]), float(line[22:30])))
                        if mgc:
                            find_elements.append(line[31:33].strip())
                        if ogc:
                            atom_info[i] = []  # Potřeba pro uložení všech vazeb jednotlivých atomů
                    if ogc:  # Příprava matic při výpočetu OGC
                        """
                        size_matrix je maximální velikost matice
                        orbital_electrons je součet valenčních elektronů jednotlivých atomu
                        count_bond_matrix je matice stupně
                        bond_matrix je matice sousednosti
                        table_elektronegativity je vektor pro tabulkové hodnoty elektronegativity
                        table_hardness je vektor pro tabulkové hodnoty tvrdosti
                        matrix_covalent_radii je vektor pro tabulkové hodnoty kovalentního poloměru
                        """
                        size_matrix, orbital_electrons = get_orbital_electrons(elements)
                        count_bond_matrix = np.zeros((size_matrix, size_matrix))
                        bond_matrix = np.zeros((size_matrix, size_matrix))
                        table_electronegativity = np.zeros((size_matrix, 1))
                        table_hardness = np.zeros((size_matrix, 1))
                        matrix_covalent_radii = np.zeros((size_matrix, 1))
                    if mgc:  # Příprava matic při výpočetu OGC
                        count_bond_matrix = np.zeros((count_atoms, count_atoms))
                        bond_matrix = np.zeros((count_atoms, count_atoms))
                    for i in range(count_bonds):  # Načtení vazeb mezi atomy
                        """
                        first_atom a second_atom je číslo atomů, mezi kterými je vazba
                        """
                        line = fh.readline()
                        first_atom, second_atom, bond = int(line[0:3]), int(line[3:6]), int(line[8:9])
                        max_bond[first_atom] = max(max_bond[first_atom], bond)
                        max_bond[second_atom] = max(max_bond[second_atom], bond)
                        if mgc:
                            """
                            index1 a index2 je vyjádření čísla atomu v maticové souřadnici
                            """
                            index1, index2 = first_atom - 1, second_atom - 1
                            bond_matrix[index1, index2] = bond_matrix[index2, index1] = bond
                            count_bond_matrix[index1, index1] += bond
                            count_bond_matrix[index2, index2] += bond
                        if ogc:
                            """
                            info_first a info_second s atom_info ukládají vazby jednolivých atomů podle pořadí
                             jak po sobě následují
                             """
                            info_first, info_second = list(atom_info[first_atom]), list(atom_info[second_atom])
                            info_first.append(bond)
                            info_second.append(bond)
                            atom_info[first_atom] = info_first
                            atom_info[second_atom] = info_second
                            orbital_electrons[0] = 0
                            position1 = position2 = 0
                            """
                            zjištění začátečních pozic atomů v matici
                            """
                            for valence in range(0, first_atom):
                                position1 += orbital_electrons[valence]
                            for valence in range(0, second_atom):
                                position2 += orbital_electrons[valence]
                            try:
                                for x in range(bond):
                                    """
                                    shift slouží k posunu na další valenční elektron, který ještě není propojen
                                    """
                                    bond_matrix[position1 + shift[first_atom]][position2 + shift[second_atom]] = \
                                        bond_matrix[position2 + shift[second_atom]][position1 + shift[first_atom]] = 1
                                    shift[first_atom] += 1
                                    shift[second_atom] += 1
                            except IndexError:
                                print("Can not prepare data from ", name)
                                break
                    try:
                        if ogc:
                            """
                            zjištění typu vazby u atomů, příprava tabulkových hodnot, přebytečných řádků
                            a sloupců v matici se provádí funkcí get_valence_state_and_prepare_table_values
                            """
                            orbital_electrons, max_bond, table_electronegativity, table_hardness,\
                            matrix_covalent_radii, valence_state, \
                            delete_rows = get_valence_state_and_prepare_table_values(atom_info, elements,
                                                                                     orbital_electrons, max_bond,
                                                                                     table_electronegativity,
                                                                                     table_hardness,
                                                                                     matrix_covalent_radii)
                            deleted = 0  # posun při vymazání řádku a sloupců
                            for row in sorted(set(delete_rows)):  # mazání v matici, nejdříve vektory
                                table_electronegativity = np.delete(table_electronegativity, row - deleted, 0)
                                table_hardness = np.delete(table_hardness, row - deleted, 0)
                                matrix_covalent_radii = np.delete(matrix_covalent_radii, row - deleted, 0)
                                for columns_or_row in [0, 1]:  # mazání rádků a sloupců v matic
                                    bond_matrix = np.delete(bond_matrix, row - deleted, columns_or_row)
                                    count_bond_matrix = np.delete(count_bond_matrix, row - deleted, columns_or_row)
                                deleted += 1
                            position = 0
                            for valence_electron in orbital_electrons:
                                """
                                zjistění zda má atom všechny atomy ve vazbě, pokud má volné elektronové páry, tak 
                                zkrátí na velikost valenční vrstvy na 4
                                doplnění vazeb v rámci atomu v matici sousednosti, místo atomu se definuje pořadím 
                                a vazby se doplňují dvěma for cykly o velikosti valenční vrstvy atomu
                                """
                                if orbital_electrons[valence_electron] >= (shift[valence_electron] + 1):
                                    while orbital_electrons[valence_electron] > 4:
                                        orbital_electrons[valence_electron] -= 1
                                for electron1 in range(orbital_electrons[valence_electron]):
                                    for electron2 in range(orbital_electrons[valence_electron]):
                                        if electron1 != electron2:
                                            bond_matrix[electron1 + position, electron2 + position] = 1
                                position += orbital_electrons[valence_electron]
                            for index in range(count_bond_matrix.shape[0]):
                                """
                                matice stupně je součet řádku v matici sousednosti 
                                """
                                count_bond_matrix[index, index] = int(sum(bond_matrix[index, :]))
                    except IndexError:
                        print("Can not prepare matrix for ", name)
                        while True:  # pokud je chyba najde se konec molekuly pro další práci v sadě
                            line = fh.readline()
                            if "$$$$" in line:
                                break
                    for number in sorted(max_bond):
                        if eem:
                            if yes_type:  # součet atomů v molekule podle prvku a maximální vazby
                                elements_count[(elements[number - 1], max_bond[number])] += 1
                            else:  # součet atomů v molekule podle prvku
                                elements_count[elements[number - 1]] += 1
                            atoms.append(Atom(number, elements[number - 1], max_bond[number], coordinate[number - 1]))
                        else:
                            if ogc:
                                """
                                elements_count se definuje při počítání ogc pro uložení velikosti valenční vrstvy
                                """
                                elements_count[elements[number - 1], number] = len(valence_state[number])
                            atoms.append(Atom(number, elements[number - 1], max_bond[number]))
                    while True:
                        line = fh.readline()
                        if "$$$$" in line:
                            if eem:
                                count_bond_matrix = bond_matrix = False  # není potřeba
                            if mgc:
                                elements_count = False  # není potřeba
                            if ogc:
                                molecules.append(Molecule(name, count_atoms, atoms, elements_count, count_bond_matrix,
                                                          bond_matrix, table_electronegativity, table_hardness,
                                                          matrix_covalent_radii))
                                """
                                ogc má vlastní způsob uložení do třídy molecules, jelikož se do ní ukládají více matic 
                                """
                            else:
                                molecules.append(Molecule(name, count_atoms, atoms, elements_count, count_bond_matrix,
                                                          bond_matrix))
                            break
        except IOError:
            print("Wrong file for molecules set! Try another file than {}".format(filename))
            sys.exit(1)

    """NAČTENÍ PARAMETRŮ PRO EEM"""
    def load_parameters(self, file_parameters):
        parameters_data, yes_type, element, kappa = [], False, "", 0.0
        try:
            with open(file_parameters, "r") as fp:
                parameters = []
                for line in fp:
                    if "Kappa=" in line:  # uložení parametru Kappa
                        kappa = float(line[line.index("Kappa") + len("Kappa=\""):-3])
                    if "<Element" in line:  # nalezení prvku pro které jsou v souboru parametry
                        element = line[line.index("Name") + len("Name=\""):-3].strip()
                    if "<Bond" in line:  # načtení parametru pro prvek o maximální vazbě
                        element_parameters = []
                        if "Type" in line:  # pokud jsou v soboru rozděleny různé hodnoty pro rozlišení typu vazeb
                            yes_type = True
                            type_bond = int(line[line.index("Type") + len("Type=\"")])
                        else:  # když není rozdělení na typy vazeb
                            type_bond = 1
                        """
                        uložení parametru A, B
                        """
                        parameter_a = int(line.index("A=") + len("A=\""))
                        parameter_b = int(line.index("B=") + len("B=\""))
                        element_parameters.append((float(line[parameter_a:parameter_a + 6]),
                                                   float(line[parameter_b:parameter_b + 6])))
                        parameters.append((element, type_bond, element_parameters))
            parameters_data.append((kappa, yes_type, parameters))
            parameter = parameters_data[0]  # jelikož jsem přehazoval moc listů pro uložení do jednoho listu
            self.parameters = parameter
            print("Load parameters for elements from {}".format(file_parameters))
        except IOError:
            print("Wrong file for parameters! Try another file than").format(file_parameters)
            sys.exit(2)

    def __str__(self):
        return str("{}".format(self.molecules))


"""FUNKCE PRO ELEKTRONEGATIVITU Z TABULEK MGC"""


def get_electronegativity_from_periodic_table(elements):
    periodic_table = {}
    for element in elements:  # list prvků, které jsou v sadě molekul
        with open("tables/Periodic Table of Elements.csv", "r", encoding="latin-1") as ftp:
            for line in ftp:
                if ("," + element + ",") in line:
                    """
                    uložení elektronegativity pro každý prvek
                    """
                    periodic_table[element] = float(line[[m.start() for m in re.finditer(r",", line)][10] + 1:[
                        m.start() for m in re.finditer(r",", line)][11]])
                    break
    return periodic_table


"""FUNKCE PRO ZÍSKANÍ TABULKOVÉ VELIKOSTI VALENČNÍ VRSTVY ATOMU"""


def get_orbital_electrons(elements):
    matrix, orbital_electrons, ionization_potential = 0, Counter(), {}
    for x, element in enumerate(elements):
        x += 1
        with open("tables/Periodic Table of Elements.csv", "r", encoding="latin-1") as ftp:
            for line in ftp:
                if ("," + element + ",") in line:  # nalezení prvku v tabulce
                    number = int(line[0:[m.start() for m in re.finditer(r",", line)][0]])  # pořadové číslov tabulce
                    if number < 3:  # když je to H nebo He
                        valence_number = int(line[[m.start() for m in re.finditer(r",", line)][-3] + 3:[
                            m.start() for m in re.finditer(r",", line)][-2]])
                        """uložení velikosti velanční vrstvy jde u těchto dvou atomů podle poslední číslice"""
                    else:
                        """
                        velikost valenční vrstvy se dělá na základě vypsané elektronové konfigurace
                        """
                        valence_number = 0
                        valence_text = line[[m.start() for m in re.finditer(r",", line)][-3] + 5:[
                            m.start() for m in re.finditer(r",", line)][-2]]
                        """
                        valence_text slouží pro uložení elektronové konfigurace
                        """

                        if len(valence_text) > 8:
                            if len(valence_text) > 10:
                                valence_text = valence_text[-8:]
                                """
                                Pokud je valence_text delší, je potřeba vzít pouze s a p konfigurace
                                """
                                if len(valence_text) > 8:
                                    valence_text = valence_text[-1]
                                    valence_number = int(valence_text[-1])
                            else:
                                """
                                valence_text je kratší než 10 znaků a delší než 8, tak lze velikost valenční vrstvy vzít 
                                z posledního čísla v zápise
                                """
                                valence_text = valence_text[-1]
                                valence_number = int(valence_text[-1])
                                break
                        maximum = int((len(valence_text)) / 4)
                        for i in range(maximum):
                            number = valence_text[3:5].strip()
                            valence_number += int(number)
                            valence_text = valence_text[3 + len(number):]
                            """
                            sčítají se orbitaly s a p a jejich zaplnění
                            """
                    orbital_electrons[x] = valence_number  # valenční velikost jednotlivých atomů
                    matrix += valence_number  # maximální velikost matic
                    break
    return matrix, orbital_electrons


"""FUNKCE PRO ELEKTRONEGATIVITU A TVRDOST Z TABULEK OGC"""


def get_electronnegativity_and_hardness(state_text, element):
    table_electronegativity, table_hardness = Counter(), Counter()
    with open("tables/electronegativity_hardness.csv", "r", encoding="latin-1") as feh:
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
                    """
                    vyhledává se v tabulce prvek a jeho zaplněná valenční vrstva a uloží se jeho elektronegativita 
                    a tvrdost pro všechny jednotlivé stavy
                    """
    return table_electronegativity, table_hardness


"""FUNKCE PRO KOVALENTNÍ POLOMĚR OGC"""


def get_covalent_radii(element, bond):
    covalent_radius = 0
    with open("tables/Covalent radii.csv", "r", encoding="latin-1") as fcr:
        for line in fcr:
            if (element + "," + str(bond)) in line:
                covalent_radius = float(line[[m.start() for m in re.finditer(r",", line)][1] + 1:])
                """
                hledání v tabulce podle symbolu prvku a jeho maximální vazby 
                """
    return covalent_radius


"""FUNKCE PRO DETEKCI VAZBY NA JEDNOTLIVÝCH ELEKTRONECH A PŘÍPRAVA MATIC"""


def get_valence_state_and_prepare_table_values(bond_info, elements, orbital_electrons, max_bond,
                                               table_electronegativity, table_hardness, matrix_covalent_radii):
    covalent_radius, valence_state, delete_rows = {}, {}, []
    orbital_electrons[0] = 0
    for number_element in bond_info:
        """ 
        uložili jsme si pro každý atom všechny jeho vazby jak jdou za sebou
        """
        position = 0
        for valence in range(0, number_element):
            position += orbital_electrons[valence]  # zjištění pozice atomu v maticích
        covalent_radius[elements[number_element - 1], max(bond_info[number_element])] = get_covalent_radii(
            elements[number_element - 1], max(bond_info[number_element]))
        """ 
        zjištění kovalentního poloměru atomu podle maximálí vazby
        """
        state_text, filled_type_state, non_binding_pair, sigma = [], 1, 0, 0
        valence_electrons = orbital_electrons[number_element]  # počet valenčních elektronů v atomu
        prepare = valence_electrons - sum(bond_info[number_element])
        long_state, index_state, text, filled = valence_electrons, 4 - (min(4, valence_electrons) - 1), "", 0
        """
        long_state je nastaven maximálně pro 4 vazbené atomy, jelikož pro více nejsou zahrnuty v použitých tabulkových
        hodnot
        """
        if prepare/2 >= 1:
            long_state = int(valence_electrons - (prepare/2))
            """ 
            prepare slouží k zjištění zda jsou všechny elektrony ve vazbě nebo se v něm nachází volný elektronový pár
            """
        if sum(bond_info[number_element]) - len(bond_info[number_element]) > 0:
            index_state = sum(bond_info[number_element]) - len(bond_info[number_element]) + 1
            """
            určení znaku pro sigma vazbu, když nemá jen sigma vazby.
            """
        electron_states = ["s", "di", "tr", "te"]  # znak pro příslušnou sigma vazbu
        divisor = 2  # dělitel
        for bond in bond_info[number_element]:  # načítání po jednotlivých vazbách v atomu
            if bond > 1:  # zapsání pi vazby
                for x in range(bond - 1):
                    state_text.append("pi")
                    filled_type_state += 1
                    filled += 1  # číslo zaplnění
                state_text.append(electron_states[-index_state])
                filled_type_state += 1
                filled += 1
                divisor = 2
                """
                zapíše se pi vazby a sigma vazba a číslo zaplnění se zvýči o počet pi vazeb a jednu sigma vazbu
                """
            else:
                state_text.append(electron_states[-index_state])
                filled += 1
                """
                zapsání sigma vazby s příslušným znakem
                """
        if ((valence_electrons - filled) / divisor) >= 1:
            for y in range((valence_electrons - filled) - (long_state - filled)):
                state_text.append(electron_states[-index_state] + "2")
                non_binding_pair += 1
                """
                po projití vazeb atomu se přiřazují volné elektronové páry, pokud se v aotmu nachází
                """
        """
        načtení tabulkových hodnot elektronegativit a tvrdosti
        """
        electronegativity, hardness = get_electronnegativity_and_hardness(state_text, elements[number_element - 1])
        valence_state[number_element] = state_text  # uložení valenčního stavu atomu
        dive = 1
        for move, state in enumerate(state_text):  # zaplnění vektorů potřebné pro finální počítání
            table_electronegativity[position + move][0] = electronegativity[
                elements[number_element - 1], state]
            table_hardness[position + move][0] = hardness[elements[number_element - 1], state]
            matrix_covalent_radii[position + move][0] = covalent_radius[elements[number_element - 1],
                                                                        max(bond_info[number_element])]
            if "2" in state:  #
                """ 
                když je v atomu volný elektronový pár, tak se zapíše pořadí posledních elektornů, jejichž množství je
                rovno počtu volných elektronových párů
                """
                if move <= 3:
                    dive += 1
                delete_rows.append(int(position + len(state_text) - 2 + dive))
    return orbital_electrons, max_bond, table_electronegativity, table_hardness, matrix_covalent_radii, valence_state,\
           delete_rows
