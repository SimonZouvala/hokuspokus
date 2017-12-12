import sys
import matplotlib.pyplot as plt
import math

from matplotlib.ticker import AutoMinorLocator, LinearLocator
from collections import Counter


class Graph:
    def __init__(self, data):
        self.first_data, self.second_data = data
        self.data1, self.data2 = load_files(self.first_data, self.second_data)
        self.ready1, self.ready2 = set_graph(self.data1, self.data2)
        graph(self.ready1, self.ready2)


def load_files(first_data, second_data):
    try:
        molecule = []
        with open(first_data, "r") as f1:
            while True:
                element_charge, line = [], f1.readline()
                if "" == line[0:1]:
                    data1 = sorted(molecule)
                    break
                name = line[:].strip()
                line = f1.readline()
                count = int(line[:].strip())
                for i in range(count):
                    line = f1.readline()
                    element = line[6:11].strip()
                    charge = float(line[11:].strip())
                    element_charge.append((element, charge))
                molecule.append((name, element_charge))
    except IOError:
        print("Wrong file for graph! Try another file than {}".format(first_data))
        sys.exit(1)
    try:
        molecule = []
        with open(second_data, "r") as f2:
            while True:
                element_charge, line = [], f2.readline()
                if "" == line[0:1]:
                    data2 = sorted(molecule)
                    break
                name = line[:].strip()
                line = f2.readline()
                count = int(line[:].strip())
                for i in range(count):
                    line = f2.readline()
                    element = line[6:11].strip()
                    charge = float(line[11:].strip())
                    element_charge.append((element, charge))
                molecule.append((name, element_charge))
    except IOError:
        print("Wrong file for graph! Try another file than {}".format(second_data))
        sys.exit(1)
    return data1, data2


def set_graph(data1, data2):
    first_data, second_data = [], []
    for name1, elements1 in data1:
        for name2, elements2 in data2:
            if name1 == name2:
                first_data.append(elements1)
                second_data.append(elements2)
                continue
    ready1 = prepare_data(first_data)
    ready2 = prepare_data(second_data)
    return ready1, ready2


def graph(ready1, ready2):
    fig, ax = plt.subplots(figsize=(12, 12))
    draw_data(ready1, ready2)
    ax.legend(loc=2)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_major_locator(LinearLocator())
    ax.xaxis.set_major_locator(LinearLocator())
    plt.tick_params(which='major', length=5, width=2)
    plt.tick_params(which='minor', length=5)
    ax.set_title("Correlation graph")
    ax.set_ylabel("Charge")
    ax.set_xlabel("Charge")
    plt.show()


def draw_data(data1, data2):
    for (element1, coordinate1), (element2, coordinate2) in zip(data1, data2):
        x, y = [], []
        for x_data in coordinate1:
            x.append(x_data)
        for y_data in coordinate2:
            y.append(y_data)
        plt.scatter(x, y, marker="o", label=(element1 + " 1. result "), alpha=0.2)
        plt.scatter(y, x, marker="2", label=(element2 + " 2. result "), alpha=0.8)
        mistake(x, y, element1)


def mistake(charges1, charges2, element):
    count, maximum, d, data_for_rmsd, suma_x, suma_y, rmsd_x, rmsd_y, covariance = 0, abs(charges1[0] - charges2[0]),\
                                                                                   0, 0, 0, 0, 0, 0, 0
    for x, y in zip(charges1, charges2):
        count += abs(x - y)
        maximum = max(maximum, abs(x - y))
        suma_x += x
        suma_y += y
    mean = count/len(charges1)
    averange_x = suma_x/(len(charges1))
    averange_y = suma_y/len(charges1)
    for x, y in zip(charges1, charges2):
        data_for_rmsd += (abs(x - y) - mean)**2
        rmsd_x += (x - averange_x)**2
        rmsd_y += (y - averange_y)**2
        covariance += (x - averange_x) * (y - averange_y)
    rmsd = math.sqrt((data_for_rmsd/len(charges1)))
    pcc = covariance/math.sqrt(rmsd_x * rmsd_y)
    print("For {} mistakes:".format(element))
    print("Mean absolute error: {}".format(mean))
    print("Maximum absolute error: {}".format(maximum))
    print("Root-mean-square deviation: {}".format(rmsd))
    print("Pearson correlation coefficient: {}". format(pcc))


def prepare_data(crude_data):
    data_ready, data = [], []
    occurrence = Counter()
    for atom in crude_data:
        for element_name, charge in atom:
            occurrence[element_name] += 1
            data.append((element_name, charge))
    data = sorted(data)
    start = 0
    for element in sorted(occurrence):
        charges = []
        for i in range(start, occurrence[element] + start):
            element_name, charge = data[i]
            charges.append(charge)
        data_ready.append((element, charges))
        start += occurrence[element]
    return data_ready
