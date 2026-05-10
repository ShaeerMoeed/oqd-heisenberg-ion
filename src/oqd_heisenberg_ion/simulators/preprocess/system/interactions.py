import shutil as sh

import numpy as np


class InteractionsFactory:
    """
    Factory for generating the required instance of the Interactions subclass. Carries a registry of Interactions subclasses

    Raises:
        Exception: if requested subclass is not found
    """

    registry = {}

    def register(cls, name, subclass):
        """
        adds the specified subclass to the registry

        Args:
            name (str): name to be used for subclass
            subclass (Type[Interactions]): Interactions subclass to be registered
        """

        cls.registry[name] = subclass

    def extract_args(cls, name, hamiltonian_term_type, **kwargs):
        """
        extracts the arguments associated with a given subclass

        Args:
            name (str): subclass name (must exist in registry)
            hamiltonian_term_type (str): Either 'zz' or 'xy'
            **kwargs (dict): key word arguments. Must contain inputs for the specific subclass

        Returns:
            (dict): contains the subclass arguments as key value pairs
        """

        arg_vals = {}
        for key, arg_dtype in cls.registry[name].args.items():
            arg_vals[key] = arg_dtype(kwargs[hamiltonian_term_type + "_" + key])

        return arg_vals

    def create(cls, hamiltonian_term_type, name, geometry, **kwargs):
        """
        creates an instance of the subclass specified

        Args:
            hamiltonian_term_type (str): either "xy" or "zz"
            name (str): name of requested subclass
            geometry (Geometry): lattice description

        Raises:
            Exception: if the requested Interactions subclass is not found in the registry

        Returns:
            (Interactions): instance of the the requested subclass
        """

        if name not in cls.registry:
            raise Exception(f"Interactions implementation not found for interaction name: {name}")
        else:
            return cls.registry[name](hamiltonian_term_type, geometry, **kwargs)


InteractionsFactory.register = classmethod(InteractionsFactory.register)
InteractionsFactory.create = classmethod(InteractionsFactory.create)
InteractionsFactory.extract_args = classmethod(InteractionsFactory.extract_args)


class Interactions:
    """
    Interactions base class. Different types of interactions implemented as subclasses
    """

    def __init__(self, hamiltonian_term_type, geometry):
        """
        constructor initializes the member variables

        Args:
            geometry (Geometry): contains contains the lattice sites and distances
        """

        self.geometry = geometry
        self.hamiltonian_term_type = hamiltonian_term_type
        self.J_ij_matrix = None
        self.J_ij_vector = None
        self.J_ij_file = None

    def get_J_ij(self):
        """
        every Interactions subclass must implement a method to populate the interaction matrix
        """

        pass

    def write_to_file(self, target_file):
        """
        every Interactions subclass must implement a method to write the interaction matrix
        """

        pass

    def update_parameters(self, parameter_dict):
        """
        updates a given parameters set with interaction properties, implemented by subclasses 

        Args:
            parameter_dict (dict): single parameter set

        Returns:
            (dict): updated parameter set
        """

        pass


class MatrixInput(Interactions):
    """
    Used when the interactions are specified via an input file containing an coupling matrix
    """

    args = {"interaction_matrix_file": str}

    def __init__(self, hamiltonian_term_type, geometry, interaction_matrix_file):
        """
        populates the J_ij vector from file

        Args:
            geometry (Geometry): contains the lattice sites and distances
            interaction_matrix_file (str): file path to the interaction matrix file
        """

        super().__init__(hamiltonian_term_type, geometry)

        self.J_ij_file = interaction_matrix_file
        geometry.initialize_tables()
        self.get_J_ij(geometry.num_bonds, geometry.sites, self.J_ij_file)

    def get_J_ij(self, num_bonds, sites, interaction_matrix_file):
        """
        populates the J_ij vector with the coupling strength for each bond

        Args:
            num_bonds (int): number of bonds in the lattice
            interaction_matrix_file (str): file path to the interactions matrix
        """

        self.J_ij_matrix = np.loadtxt(interaction_matrix_file, delimiter=",", skiprows=1)
        self.J_ij_vector = np.zeros(num_bonds)
        for b in range(num_bonds):
            self.J_ij_vector[b] = self.J_ij_matrix[sites[b,0], sites[b,1]]

    def write_to_file(self, target_file):
        """
        copies the input J_ij matrix to specified directory

        Args:
            target_file (str): path to target file
        """

        sh.copyfile(self.J_ij_file, target_file)

    def update_parameters(self, parameter_dict):
        """
        updates a given parameters set with interaction properties

        Args:
            parameter_dict (dict): single parameter set

        Returns:
            (dict): updated parameter set
        """

        parameter_dict[self.hamiltonian_term_type + "_interaction_type"] = "matrix_input"
        parameter_dict[self.hamiltonian_term_type + "_interaction_matrix_file"] = self.J_ij_file

        return parameter_dict


class VectorInput(Interactions):
    """
    Used when the interactions are specified via an input file containing an interaction vector indexed by bond indices
    """

    args = {"interaction_vector_file": str}

    def __init__(self, hamiltonian_term_type, geometry, interaction_vector_file):
        """
        populates the J_ij vector from file

        Args:
            geometry (Geometry): contains the lattice sites and distances
            interaction_vector_file (str): file path to the interaction vector file
        """

        super().__init__(hamiltonian_term_type, geometry)

        self.J_ij_file = interaction_vector_file
        geometry.initialize_tables()
        self.get_J_ij(geometry.N, geometry.num_bonds, geometry.sites, self.J_ij_file)

    def get_J_ij(self, N, num_bonds, sites, interaction_vector_file):
        """
        populates the J_ij vector with the coupling strength for each bond

        Args:
            num_bonds (int): number of bonds in the lattice
            interaction_vector_file (str): file path to the interactions vector
        """

        bonds, J_ij_vector_from_file = np.loadtxt(interaction_vector_file, delimiter=",", skiprows=1, unpack=True)
        if len(J_ij_vector_from_file != num_bonds):
            raise Exception("The interaction vector provided via file is not consistent with geometry")
        self.J_ij_matrix = np.zeros((N, N))
        for b in range(num_bonds):
            self.J_ij_matrix[sites[b,0], sites[b,1]] = self.J_ij_vector[b]

    def write_to_file(self, target_file):
        """
        writes the input J_ij matrix to specified directory

        Args:
            target_file (str): path to target file
        """

        np.savetxt(
            target_file,
            self.J_ij_matrix,
            delimiter=",",
            header="J_ij_matrix, interactions=vector_input",
        )

    def update_parameters(self, parameter_dict):
        """
        updates a given parameters set with interaction properties 

        Args:
            parameter_dict (dict): single parameter set

        Returns:
            (dict): updated parameter set
        """

        parameter_dict[self.hamiltonian_term_type + "_interaction_type"] = "vector_input"
        parameter_dict[self.hamiltonian_term_type + "_interaction_matrix_file"] = self.J_ij_file

        return parameter_dict


class PowerLaw(Interactions):
    """
    Used to generate the coupling using power law J_ij = 1/r_{ij}^alpha
    """

    args = {"alpha": float}

    def __init__(self, hamiltonian_term_type, geometry, alpha):
        """
        constructs the coupling matrix

        Args:
            geometry (Geometry): contains the lattice sites and distances
            alpha (float): power law interaction strength exponent
        """

        super().__init__(hamiltonian_term_type, geometry)

        self.alpha = alpha
        self.geometry.initialize_tables()
        self.geometry.build()
        self.get_J_ij(self.geometry.N, self.geometry.num_bonds, self.geometry.distances, self.geometry.sites, self.alpha)

    def get_J_ij(self, N, num_bonds, distances, sites, alpha):
        """
        populates the J_ij matrix

        Args:
            N (int): number of sites
            num_bonds (int): number of bonds
            distances (numpy.ndarray[float]): num_bonds x 1 array containing the distances between all pairs of interacting sites
            sites (numpy.ndarray[float]): numbonds x 2 array containing the site indices corresponding to all bonds
            alpha (float): interaction strength exponent
        """

        self.J_ij_vector = np.zeros(num_bonds)
        self.J_ij_matrix = np.zeros((N, N))
        for b in range(num_bonds):
            r_pow_alpha = (distances[b]) ** alpha
            self.J_ij_matrix[sites[b,0], sites[b,1]] = 1.0 / r_pow_alpha
            self.J_ij_matrix[sites[b,1], sites[b,0]] = 1.0 / r_pow_alpha
            self.J_ij_vector[b] = 1.0 / r_pow_alpha

    def write_to_file(self, target_file):
        """
        writes the J_ij matrix to file

        Args:
            target_file (str): path to target file
        """

        np.savetxt(
            target_file,
            self.J_ij_matrix,
            delimiter=",",
            header=f"J_ij_matrix, interactions=power_law,alpha={self.alpha}",
        )

    def update_parameters(self, parameter_dict):
        """
        updates a given parameters set with interaction properties

        Args:
            parameter_dict (dict): single parameter set

        Returns:
            (dict): updated parameter set
        """

        parameter_dict[self.hamiltonian_term_type + "_interaction_type"] = "power_law"
        parameter_dict[self.hamiltonian_term_type + "_alpha"] = self.alpha

        return parameter_dict


class Constant(Interactions):
    """
    Constant interactions
    """

    args = {}

    def __init__(self, hamiltonian_term_type, geometry):
        """
        initializes constant interactions

        Args:
            geometry (Geometry): contains the lattice sites and distances
        """

        super().__init__(hamiltonian_term_type, geometry)

        self.geometry.initialize_tables()
        self.geometry.build()
        self.get_J_ij(self.geometry.N, self.geometry.num_bonds, self.geometry.sites)

    def get_J_ij(self, N, num_bonds, sites):
        """
        populates the J_ij matrix with J_ij=1 for all bonds

        Args:
            N (int): number of sites
            num_bonds (int): number of bonds
            sites (numpy.ndarray[float]): numbonds x 2 array containing the site indices corresponding to all bonds
        """

        N = self.geometry.N
        num_bonds = self.geometry.num_bonds
        sites = self.geometry.sites
        self.J_ij_vector = np.zeros(num_bonds)
        self.J_ij_matrix = np.zeros((N, N))
        for b in range(num_bonds):
            self.J_ij_matrix[sites[b, 0], sites[b, 1]] = 1.0
            self.J_ij_matrix[sites[b, 1], sites[b, 0]] = 1.0
            self.J_ij_vector[b] = 1.0

    def write_to_file(self, target_file):
        """
        writes the J_ij matrix to file

        Args:
            target_file (str): path to target file
        """

        np.savetxt(
            target_file,
            self.J_ij_matrix,
            delimiter=",",
            header="J_ij_matrix, interactions=constant",
        )

    def update_parameters(self, parameter_dict):
        """
        updates a given parameters set with interaction properties

        Args:
            parameter_dict (dict): single parameter set

        Returns:
            (dict): updated parameter set
        """

        parameter_dict[self.hamiltonian_term_type + "_interaction_type"] = "constant"

        return parameter_dict

InteractionsFactory.register("power_law", PowerLaw)
InteractionsFactory.register("matrix_input", MatrixInput)
InteractionsFactory.register("vector_input", VectorInput)
InteractionsFactory.register("constant", Constant)
