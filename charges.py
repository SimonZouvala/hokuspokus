class Charges:
    def __init__(self, filename, file_parameters, data):
        data_for_graph = []
        self.new_file = "result/data_from_" + filename.strip("data/") + "_with_parameters_"\
                        + file_parameters.strip("data/")
        with open(self.new_file, "w") as f:
            for name, count, atoms, charges, parameters in data:
                try:
                    print("{}\n{}".format(name, int(count)), file=f)
                    for i, atom in enumerate(atoms):
                        element, number, bond = atom
                        print("{0:6d}  {1:>2}{2} {3: f}".format(number, element, bond, charges[i]), file=f)
                    data_for_graph.append((name, count, atoms, charges, parameters))
                except ValueError:
                    continue
        print("Now you can find charge for each element in file {}".format(self.new_file))

    def get_name_file_with_result(self):
        return self.new_file
