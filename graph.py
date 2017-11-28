import sys
import argparse
import math
import numpy as np
import matplotlib.pyplot as plt
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
        draw_data(self.ready1, "o", 1)
        draw_data(self.ready2, ",", 2)
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


def draw_data(data_for_draw, marker, number):
    for element, coordinate in sorted(data_for_draw):
        x, y,  = [], []
        for y_data in coordinate:
            y.append(y_data)
        plt.scatter(y, y, marker=marker, label=(element + " from " + str(number) + ". result "), alpha=0.5)


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
        prepare_data = []
        for i in range(start, occurrence[element] + start):
            element_name, charge = data[i]
            prepare_data.append(charge)
        data_ready.append((element, prepare_data))
        start += occurrence[element]
    return data_ready


draw = Graph("result/data_from_set01.sdf_with_parameters_ElemBond.par",
             "result/data_from_set01.sdf_with_parameters_Element.par")
draw.load_files()
draw.set_graph()
draw.graph()
