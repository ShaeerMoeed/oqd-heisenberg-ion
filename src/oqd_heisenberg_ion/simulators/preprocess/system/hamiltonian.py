import numpy as np


class HamiltonianParameters:
    """
    Base class for Hamiltonian parameters
    """

    args = {}

    def __init__(self):
        """
        initializes member variables
        """

        self.hamiltonian_name = None
        self.hamiltonian_type = None
        self.Delta = None
        self.h = None
        self.J = None
        self.B = None
        self.xy_coeff_array = None
        self.zz_coeff_array = None
        self.z_field_array = None


    def update_parameters(self, parameter_dict):
        """
        updates the provided dictionary with hamiltonian parameters

        Args:
            parameter_dict (dict): dict to be updated

        Returns:
            (dict): updated dict containing the parameter set
        """

        parameter_dict["hamiltonian_name"] = self.hamiltonian_name
        parameter_dict["hamiltonian_type"] = self.hamiltonian_type
        parameter_dict["Delta"] = self.Delta
        parameter_dict["h"] = self.h
        parameter_dict["J"] = self.J
        parameter_dict["B"] = self.B

        return parameter_dict
    
    def construct_coeff_arrays(self, xy_interactions, zz_interactions, field):
        """
        initializes the arrays containing the coefficients of the (S^X_iS^X_j + S^Y_iS^Y_j) and S^Z_i S^Z_j terms in the Hamiltonian.
        Base implementation can be used when Delta is bond-independent.  

        Args:
            geometry (Geometry): specifies the lattice sites and distances
            xy_interactions (Interactions): contains the J_ij interactions for each bond in the lattice for the XY term
            zz_interactions (Interactions): contains the J_ij interactions for each bond in the lattice for the ZZ term
            bond_field_prefactor (float): pre-factor for field strength so that h(bond) = h * bond_field_prefactor
        """

        self.xy_coeff_array = xy_interactions.J_ij_vector
        self.zz_coeff_array = self.Delta * zz_interactions.J_ij_vector
        self.z_field_array = (self.h/self.J) * field.H_b_vector

    def validate_coefficient_arrays(self):
        """
        Each subclass should implement its own validation logic for the coefficient arrays wherever required 
        """
        pass


class HamiltonianFactory:
    """
    Factory for generating the required instance of the HamiltonianParameters subclass. Carries a registry of HamiltonianParameters subclasses

    Raises:
        Exception: if requested subclass is not found
    """

    registry = {}

    def register(cls, name, subclass):
        """
        adds the specified subclass to the registry

        Args:
            name (str): name to be used for subclass
            subclass (Type[HamiltonianParameters]): HamiltonianParameters subclass to be registered
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

    def create(cls, name, **kwargs):
        """
        creates an instance of the subclass specified

        Args:
            name (str): name of requested subclass

        Raises:
            Exception: if the requested HamiltonianParameters subclass is not found in the registry

        Returns:
            (HamiltonianParameters): instance of the the requested subclass
        """

        if name not in cls.registry:
            raise Exception(f"HamiltonianParameters implementation not found for name: {name}")
        else:
            return cls.registry[name](**kwargs)


# Make register and create classmethods so subclasses can be added to the registry
# and instantiated agnostically
HamiltonianFactory.register = classmethod(HamiltonianFactory.register)
HamiltonianFactory.create = classmethod(HamiltonianFactory.create)
HamiltonianFactory.extract_args = classmethod(HamiltonianFactory.extract_args)


class FMHeisenbergAFMZ(HamiltonianParameters):
    """
    ferromagnetic XXZ model with Delta = -1
    """

    args = {"J": float}

    def __init__(self, J):
        """
        sets the corresponding member variables

        Args:
            J (float): energy scale, must be positive for ferromagnetic interactions

        Raises:
            Exception: if J <= 0
        """

        super().__init__()

        self.hamiltonian_name = "fm_heisenberg_afm_Z"
        self.hamiltonian_type = -1
        self.Delta = -1.0
        self.h = 0.0
        self.J = J
        self.B = 0.0

        if self.J <= 0:
            raise Exception("J must be positive for ferromagnetic interactions")
        
    def validate_coefficient_arrays(self):

        if np.any(np.abs(self.z_field_array) > 1e-15):
            raise Exception("Z field coefficients should be identically zero for the FMHeisenbergAFMZ model")
        
        if np.any(np.abs(self.xy_coeff_array + self.zz_coeff_array) > 1e-15):
            raise Exception("ZZ coefficients should be equal to -XX coefficients for the FMHeisenbergAFMZ model")


class XY(HamiltonianParameters):
    """
    XY model
    """

    args = {"J": float}

    def __init__(self, J):
        """
        sets the corresponding member variables

        Args:
            J (float): energy scale
        """

        super().__init__()

        self.hamiltonian_name = "XY"
        self.hamiltonian_type = 0
        self.Delta = 0.0
        self.h = 0.0
        self.J = J
        self.B = 0.0

    def validate_coefficient_arrays(self):

        if np.any(np.abs(self.z_field_array) > 1e-15):
            raise Exception("Z field coefficients should be identically zero for the XY model")
        
        if np.any(np.abs(self.zz_coeff_array) > 1e-15):
            raise Exception("ZZ coefficients should be identically zero for the XY model")


class FMHeisenbergFMZ(HamiltonianParameters):
    """
    ferromagnetic Heisenberg model

    Raises:
        Exception: if input energy scale J <= 0
    """

    args = {"J": float}

    def __init__(self, J):
        """
        sets the corresponding member variables

        Args:
            J (float): energy scale, must be positive for ferromagnetic interactions

        Raises:
            Exception: if J <= 0
        """

        super().__init__()

        self.hamiltonian_name = "fm_heisenberg_fm_Z"
        self.hamiltonian_type = 1
        self.Delta = 1.0
        self.h = 0.0
        self.J = J
        self.B = 0.0

        if self.J <= 0:
            raise Exception("J must be positive for ferromagnetic interactions")
        
    def validate_coefficient_arrays(self):

        if np.any(np.abs(self.z_field_array) > 1e-15):
            raise Exception("Z field coefficients should be identically zero for the FMHeisenbergFMZ model")
        
        if np.any(np.abs(self.xy_coeff_array - self.zz_coeff_array) > 1e-15):
            raise Exception("ZZ coefficients should be equal to XX coefficients for the FMHeisenbergFMZ model")


class XXZ(HamiltonianParameters):
    """
    XXZ model
    """

    args = {"Delta": float, "J": float}

    def __init__(self, Delta, J):
        """
        sets the corresponding member variables

        Args:
            Delta (float): coefficient of the Z_i Z_j term in the Hamiltonian
            J (float): energy scale
        """

        super().__init__()

        self.hamiltonian_name = "XXZ"
        self.Delta = Delta
        self.hamiltonian_type = 2
        self.h = 0.0
        self.J = J
        self.B = 0.0

    def validate_coefficient_arrays(self):

        if np.any(np.abs(self.z_field_array) > 1e-15):
            raise Exception("Z field coefficients should be identically zero for the XXZ model")
        
        if np.any(np.abs(self.xy_coeff_array - (1.0/self.Delta) * self.zz_coeff_array) > 1e-15):
            raise Exception("ZZ coefficients should be equal to Delta x XX coefficients for the XXZ model")


class XXZh(HamiltonianParameters):
    """
    XXZ model with a longitudinal field
    """

    args = {"Delta": float, "h": float, "J": float}

    def __init__(self, Delta, h, J):
        """
        sets the corresponding member variables

        Args:
            Delta (float): coefficient of the Z_i Z_j term in the Hamiltonian
            h (float): longitudinal field strength, must be positive for QMC
            J (float): energy scale
        """

        super().__init__()

        self.hamiltonian_name = "XXZh"
        self.Delta = Delta
        self.hamiltonian_type = 3
        self.h = h
        self.J = J
        self.B = 0.0

    def validate_coefficient_arrays(self):

        if not np.allclose(self.z_field_array, self.z_field_array[0]):
            raise Exception("Z field coefficients should be identical for the XXZh model")
        
        if np.any(np.abs(self.xy_coeff_array - (1.0/self.Delta) * self.zz_coeff_array) > 1e-15):
            raise Exception("ZZ coefficients should be equal to Delta x XX coefficients for the XXZh model")


class AFMHeisenbergFMZ(HamiltonianParameters):
    """
    anti-ferromagnetic Heisenberg model

    Raises:
        Exception: if input energy scale J >= 0
    """

    args = {"J": float}

    def __init__(self, J):
        """
        sets the corresponding member variables

        Args:
            J (float): energy scale, must be negative for anti-ferromagnetic interactions

        Raises:
            Exception: if J >= 0
        """

        super().__init__()

        self.hamiltonian_name = "afm_heisenberg_fm_Z"
        self.hamiltonian_type = 4
        self.Delta = 1.0
        self.h = 0.0
        self.J = J
        self.B = 0.0

        if self.J >= 0:
            raise Exception("J must be negative for anti-ferromagnetic interactions")
        
    def validate_coefficient_arrays(self):

        if np.any(np.abs(self.z_field_array) > 1e-15):
            raise Exception("Z field coefficients should be identically zero for the AFMHeisenbergFMZ model")
        
        if np.any(np.abs(self.xy_coeff_array - self.zz_coeff_array) > 1e-15):
            raise Exception("ZZ coefficients should be equal to XX coefficients for the AFMHeisenbergFMZ model")


class XXZhB(HamiltonianParameters):
    """
    XXZ model with a longitudinal field and a transverse field. Only usable with ED
    """

    args = {"Delta": float, "h": float, "B": float, "J": float}

    def __init__(self, Delta, h, B, J):
        """
        sets the corresponding member variables

        Args:
            Delta (float): coefficient of the Z_i Z_j term in the Hamiltonian
            h (float): longitudinal field strength
            B (float): transverse field strength
            J (float): energy scale
        """

        super().__init__()

        self.hamiltonian_name = "XXZh"
        self.Delta = Delta
        self.hamiltonian_type = 5
        self.h = h
        self.J = J
        self.B = B

    def validate_coefficient_arrays(self):

        if not np.anyclose(self.z_field_array, self.z_field_array[0]):
            raise Exception("Z field coefficients should be identical for the XXZhB model")
        
        if np.any(np.abs(self.xy_coeff_array - self.zz_coeff_array) > 1e-15):
            raise Exception("ZZ coefficients should be equal to Delta x XX coefficients for the XXZhB model")
        

class InhomogenousXXZh(HamiltonianParameters):

    """
    Inhomogenous XXZ model with a longitudinal field. Allows independent ZZ and XY coefficients and site-dependent fields
    TODO: Add support for long range QMC and ED for this hamiltonian type
    """

    args = {"Delta": float, "h": float, "J": float}

    def __init__(self, Delta, h, J):
        """
        sets the corresponding member variables

        Args:
            Delta (float): coefficient of the Z_i Z_j term in the Hamiltonian
            h (float): longitudinal field strength
            J (float): energy scale
        """

        super().__init__()

        self.hamiltonian_name = "inhomogenous_XXZh"
        self.Delta = Delta
        self.hamiltonian_type = 6
        self.h = h
        self.J = J


# Register all implemented Hamiltonian parameter types in the base HamiltonianParameters class
HamiltonianFactory.register("XY", XY)
HamiltonianFactory.register("XXZ", XXZ)
HamiltonianFactory.register("XXZh", XXZh)
HamiltonianFactory.register("XXZhB", XXZhB)
HamiltonianFactory.register("fm_heisenberg_fm_Z", FMHeisenbergFMZ)
HamiltonianFactory.register("fm_heisenberg_afm_Z", FMHeisenbergAFMZ)
HamiltonianFactory.register("afm_heisenberg_fm_Z", AFMHeisenbergFMZ)
HamiltonianFactory.register("inhomogenous_XXZh", InhomogenousXXZh)
