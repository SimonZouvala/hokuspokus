import sys  # pro ukončení práce při chybném souboru
import matplotlib.pyplot as plt  # knihovna Matplotlib
import math  # knihovna pro použití matematických funkcí
import numpy as np  # knihovna NumPy

from matplotlib.ticker import AutoMinorLocator, LinearLocator  # knihovna na vygenerování os x a y
from collections import defaultdict  # knihovna pro použití funkce Counter()


class Graph:
    def __init__(self, data):
        self.first_data, self.second_data = data
        self.data1, self.data2 = load_files(self.first_data), load_files(self.second_data)  # načtení dat
        self.ready1, self.ready2 = compare_data(self.data1, self.data2)  # porovnání molekul z obou sad
        graph(self.ready1, self.ready2, data)  # vykreslení grafu a spočítání statistiky


"""NAČTENÍ DAT"""


def load_files(file):
    try:
        molecule = []
        with open(file, "r") as f1:
            while True:
                """načtení souboru s náboji (načítá uložené soubory z výpočtů ogcm.py, mgcm.py nebo eem.py)"""
                element_charge, line = [], f1.readline()
                if "" == line[0:1]:
                    data = sorted(molecule)
                    break
                name = line[:].strip()
                line = f1.readline()
                count = int(line[:].strip())
                """zjištění názvu, počtu atomů, jednotlivých atomů s maximální vazbou a jejich nábojů """
                for i in range(count):
                    line = f1.readline()
                    element = line[6:11].strip()
                    charge = float(line[11:].strip())
                    element_charge.append((element, charge))
                molecule.append((name, element_charge))
    except IOError:
        print("Wrong file for graph! Try another file than {}".format(file))
        sys.exit(1)
    return data


"""POROVNÁNÍ MOLEKUL Z OBOU SAD"""


def compare_data(data1, data2):
    first_data, second_data = [], []
    """
    Zjištění zda sady obsahují stejné molekuly a ty se uloží do porovnání. Pokud je molekula, která je jen v 
    jedné sadě, tak je vynechána
    """
    for name1, elements1 in data1:
        for name2, elements2 in data2:
            if name1 == name2:
                first_data.append(elements1)
                second_data.append(elements2)
                continue
    ready1 = prepare_data(first_data)  # příprava dat ze sad
    ready2 = prepare_data(second_data)
    return ready1, ready2


"""VYKRESLENÍ GRAFU A SPOČÍTÁNÍ STATISTIKY"""


def graph(ready1, ready2, filenames):
    fig, ax = plt.subplots(figsize=(12, 12))  # velikost vykresleného grafu
    minimum, maximum = draw_data(ready1, ready2)  # vykreslení bodů v grafu a zjištění extrémů
    ax.legend(loc=2)  # kde se zobrazí legenda
    """vygenerování os x a y"""
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_major_locator(LinearLocator())
    ax.xaxis.set_major_locator(LinearLocator())
    """nastavení značek osách"""
    plt.tick_params(which='major', length=5, width=2)
    plt.tick_params(which='minor', length=5)
    """nastavení stejných os"""
    plt.xlim(minimum - 0.05, maximum + 0.05)
    plt.ylim(minimum - 0.05, maximum + 0.05)
    """nastavneí popisu os"""
    ax.set_title("Correlation graph")
    ax.set_ylabel("Charge from {}".format(filenames[1]))
    ax.set_xlabel("Charge from {}".format(filenames[0]))
    plt.show()


"""VYKRESLENÍ BODŮ V GRAFU A ZJIŠTĚNÍ EXTRÉMŮ"""


def draw_data(element_data1, element_data2):
    statistics, minimum, maximum, yminimum, ymaximum = [], 99999, -99999, 99999, -999999
    pcc_graph = get_pcc_for_all_graph(element_data1, element_data2)  # zjištění celkové korelace
    for element1, element2 in zip(element_data1, element_data2):
        """jednotlivé atomy a všechny jejich náboje se používají pro zjištění extrémů """
        yminimum = min(min(element_data2[element2]), yminimum)  # extremy pro druhou sadu, pokud potřebujeme
        ymaximum = max(max(element_data2[element2]), ymaximum)
        """extremy pro obě sady"""
        minimum = min(minimum, min(element_data1[element1]), min(element_data2[element2]))
        maximum = max(maximum, max(element_data1[element1]), max(element_data2[element2]))
        """vykreslení atomu a jeho náboje"""
        plt.scatter(element_data1[element1], element_data2[element2], marker="o", label=element1, alpha=0.9)
        statistics.append(get_statistics(element_data1[element1], element_data2[element2], element1))  # statistika
    print_statistics(statistics, pcc_graph)  # vypsat statistiku
    return minimum, maximum


"""VYPSAT STATISTIKU"""


def print_statistics(statistics, pcc_graph):
    """vypíše statistiku v terminálu."""
    print("Element          MAE      ABSMAX     RMSD        PCC")
    for element, mae, maximum, rmsd, pcc in statistics:
        print("{:<10} {: >9.3} {: >11.3} {: >8.3} {: >10.3}".format(element, mae, maximum, rmsd, pcc))
    print("r = {}". format(pcc_graph))


"""STATISTIKA"""


def get_statistics(charges1, charges2, element):
    maximum, data_for_rmsd, sum_x, sum_y = abs(charges1[0] - charges2[0]), 0, 0, 0
    charges = charges1 + charges2
    value = 0
    data_x, data_y = [], []
    """procházíme hodnoty nábojů z obou sad pro zjištění průměrné hodnoty"""
    for x, y in zip(charges1, charges2):
        maximum = max(maximum, abs(x - y))  # výpočet ABSMAX
        sum_x += x
        sum_y += y
    mean = (sum_x+sum_y)/(len(charges1)+len(charges2))  # průměr
    """výpočet MAE"""
    for charge in charges:
        value += abs(charge - mean)
    mae = value/len(charges)
    """výpočet RMSD"""
    for x, y in zip(charges1, charges2):
        data_for_rmsd += (x - y)**2
        data_x.append(x)
        data_y.append(y)
    rmsd = math.sqrt((data_for_rmsd/len(charges1)))
    """výpočet PCC"""
    try:
        pcc = np.corrcoef(data_x, data_y)[1, 0]  # výpočet pomocí knihovny NumPy
    except ArithmeticError:
        pcc = 'nan'
        pass
    return [element, mae, maximum, rmsd, pcc]


"""PŘÍPRAVA DAT ZE SAD"""


def prepare_data(crude_data):
    element_data, elements = defaultdict(list), []
    """uložení všech nábojů v molekule, které přísluší danému atomu"""
    for molecule in crude_data:
        elements = [element for element, _ in molecule]
        for e in set(elements):
            element_data[e].extend(list([charge for element, charge in molecule if element == e]))
    return element_data


"""ZJIŠTĚNÍ CELKOVÉ KORELACE"""


def get_pcc_for_all_graph(charges1, charges2):
    data_x, data_y = [], []
    """potřeba sloučit všechny náboje do dvou stejně velkých polí"""
    for i, j in zip(charges1, charges2):
        for x, y in zip(charges1[i], charges2[j]):
            data_x.append(x)
            data_y.append(y)
    pearr = np.corrcoef(data_x, data_y)[1, 0]  # výpočet pomocí knihovny NumPy
    return pearr
