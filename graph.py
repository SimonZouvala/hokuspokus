import sys
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator, LinearLocator
from collections import Counter


class Graph:
    def __init__(self, first_data, second_data):
        self.first_data = first_data
        self.second_data = second_data

    def load_files(self):
        try:
            molecule = []
            with open(self.first_data, "r") as f1:
                while True:
                    element_charge, line = [], f1.readline()
                    if "" == line[0:1]:
                        self.data1 = sorted(molecule)
                        break
                    name = line[:].strip()
                    line = f1.readline()
                    count = int(line[:].strip())
                    for i in range(count):
                        line = f1.readline()
                        element = line[6:11].strip()
                        charge = float(line[11:].strip())
                        element_charge.append((element, charge))
                    molecule.append((name, sorted(element_charge)))
        except IOError:
            print("Wrong file for graph! Try another file than {}".format(self.first_data))
            sys.exit(1)
        try:
            molecule = []
            with open(self.second_data, "r") as f2:
                while True:
                    element_charge, line = [], f2.readline()
                    if "" == line[0:1]:
                        self.data2 = sorted(molecule)
                        break
                    name = line[:].strip()
                    line = f2.readline()
                    count = int(line[:].strip())
                    for i in range(count):
                        line = f2.readline()
                        element = line[6:11].strip()
                        charge = float(line[11:].strip())
                        element_charge.append((element, charge))
                    molecule.append((name, sorted(element_charge)))
        except IOError:
            print("Wrong file for graph! Try another file than {}".format(self.second_data))
            sys.exit(1)

    def set_graph(self):
        first_data, second_data = [], []
        for name1, elements1 in self.data1:
            for name2, elements2 in self.data2:
                if name1 == name2:
                    first_data.append(elements1)
                    second_data.append(elements2)
                    continue
        self.ready1 = prepare_data(first_data)
        self.ready2 = prepare_data(second_data)

    def graph(self):
        fig, ax = plt.subplots(figsize=(12, 12))
        draw_data(sorted(self.ready1), sorted(self.ready2))
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
        plt.scatter(x, y, marker="o", label=(element1 + " 1. result "), alpha=0.5)
        plt.scatter(y, x, marker=",", label=(element2 + " 2. result "), alpha=0.5)


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
