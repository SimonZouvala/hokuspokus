class Atom:
    def __init__(self, number, element_symbol, bond, coordinate):
        self.element_symbol = element_symbol
        self.coordinate = coordinate
        self.bond = bond
        self.number = number

    def __str__(self):
        return str("{} {} {} {}".format(self.number, self.element_symbol, self.bond, self.coordinate))


class Molecule:

    def __init__(self, name, count_atoms, atoms, elements_count):
        self.name = name
        self.count_atoms = count_atoms
        self.atoms = atoms
        self.elements_count = elements_count

    def __str__(self):
        return str("Molecule name: {}\n{}".format(self.name, self.atoms))


class MoleculesSet:
    def __init__(self, molecules):
        self.molecules = molecules

    def __str__(self):
        return str("{}".format(self.molecules))


class Parameters:
    def __init__(self, parameters):
        self.parameters = parameters
