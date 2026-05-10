import shutil as sh

import numpy as np


class FieldFactory:
    """
    Factory for generating the required instance of the (longitudinal) Field subclass. Carries a registry of Field subclasses

    Raises:
        Exception: if requested subclass is not found
    """

    registry = {}

    def register(cls, name, subclass):
        """
        adds the specified subclass to the registry

        Args:
            name (str): name to be used for subclass
            subclass (Type[Field]): Field subclass to be registered
        """

        cls.registry[name] = subclass

    def extract_args(cls, name, **kwargs):
        """
        extracts the arguments associated with a given subclass

        Args:
            name (str): subclass name (must exist in registry)
            **kwargs (dict): key word arguments. Must contain inputs for the specific subclass

        Returns:
            (dict): contains the subclass arguments as key value pairs
        """

        arg_vals = {}
        for key, arg_dtype in cls.registry[name].args.items():
            arg_vals[key] = arg_dtype(kwargs[key])

        return arg_vals

    def create(cls, name, geometry, **kwargs):
        """
        creates an instance of the subclass specified

        Args:
            name (str): name of requested subclass

        Raises:
            Exception: if the requested Field subclass is not found in the registry

        Returns:
            (Field): instance of the the requested subclass
        """

        if name not in cls.registry:
            raise Exception(f"Field implementation not found for field name: {name}")
        else:
            return cls.registry[name](geometry, **kwargs)


FieldFactory.register = classmethod(FieldFactory.register)
FieldFactory.create = classmethod(FieldFactory.create)
FieldFactory.extract_args = classmethod(FieldFactory.extract_args)


class Field:
    """
    Field base class. Different types of (longitudinal) fields implemented as subclasses. Specifies a vector H_b such that
    the field contribution of each bond h_b = (h/J) * H_b_vector
    """

    def __init__(self, geometry):
        """
        constructor initializes the member variables

        Args:
            geometry (Geometry): contains contains the lattice sites and distances
        """

        self.geometry = geometry
        self.H_b_vector = None
        self.H_b_file = None

    def get_H_b(self):
        """
        every Field subclass must implement a method to populate the field vector
        """

        pass

    def write_H_b_file(self, target_file):
        """
        every Field subclass must implement a method to write the field file
        """


class FileInput(Field):
    """
    Used when the field is specified via an input file
    """

    args = {"z_field_vector_file": str}

    def __init__(self, geometry, H_b_file):
        """
        populates the H_b vector from file

        Args:
            geometry (Geometry): contains the lattice sites and distances
            H_b_file (str): file path to the field file
        """

        super().__init__(geometry)

        self.H_b_file = H_b_file
        geometry.initialize_tables()
        self.get_H_b(geometry.num_bonds, self.H_b_file)

    def get_H_b(self, num_bonds, H_b_file):
        """
        populates the H_b vector with the field strength for each bond

        Args:
            num_bonds (int): number of bonds in the lattice
            H_b_file (str): file path to the file containing the field contributions for each bond
        """

        bonds_from_file, H_b_vector_from_file = np.loadtxt(H_b_file, delimiter=",", skiprows=1, unpack=True)
        if (len(H_b_vector_from_file) != num_bonds):
            raise Exception("File input for H_b vector is inconsistent with geometry")
        else:
            self.H_b_vector = H_b_vector_from_file

    def write_to_file(self, target_file):
        """
        copies the input H_b file to specified directory

        Args:
            target_file (str): path to target file
        """

        sh.copyfile(self.H_b_file, target_file)


class Constant(Field):
    """
    constant field for all sites
    """

    args = {}

    def __init__(self, geometry):
        """
        constructs the H_b vector for a constant field strength

        Args:
            geometry (Geometry): contains the lattice sites and distances
        """

        super().__init__(geometry)

        geometry.initialize_tables()
        self.get_H_b(geometry.num_bonds, geometry.num_neighbors_per_site)

    def get_H_b(self, num_bonds, num_neighbors):
        """
        populates the H_b vector with a constant field strength of 1/(num_neighbors) for each bond. 
        The bond field contribution should then be h/J * H_b_vector

        Args:
            num_bonds (int): number of bonds in the lattice
        """
        
        self.H_b_vector = (1.0/num_neighbors) * np.ones(num_bonds)

    def write_to_file(self, target_file):
        """
        writes the H_b vector to a file

        Args:
            target_file (str): path to target file
        """

        indices = np.arange(self.geometry.num_bonds)
        result = np.column_stack((indices, self.H_b_vector))

        np.savetxt(
            target_file,
            result,
            delimiter=",",
            header="H_b_vector, interactions=constant \n bond, H_b",
        )


FieldFactory.register("file_input", FileInput)
FieldFactory.register("constant", Constant)
