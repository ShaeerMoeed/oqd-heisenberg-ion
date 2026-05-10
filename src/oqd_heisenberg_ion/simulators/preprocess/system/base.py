from .field import FieldFactory
from .geometry import GeometryFactory
from .hamiltonian import HamiltonianFactory
from .interactions import InteractionsFactory


class System:
    """
    Defines the model. This includes the Hamiltonian parameters, interaction strengths and the lattice geometry
    """

    def __init__(self, **kwargs):
        """
        Defines the system from key word arguments

        Args:
            **kwargs(dict): key word arguments required to specify the system
        """

        self.model_name = kwargs["hamiltonian_name"]
        hamiltonian_args = HamiltonianFactory.extract_args(self.model_name, **kwargs)
        self.hamiltonian_parameters = HamiltonianFactory.create(self.model_name, **hamiltonian_args)

        self.geometry_name = (
            kwargs["interaction_range"]
            + "_"
            + kwargs["boundary"]
            + "_"
            + kwargs["spatial_dimension"]
            + "_"
            + kwargs["lattice_type"]
        )
        geometry_args = GeometryFactory.extract_args(self.geometry_name, **kwargs)
        self.geometry = GeometryFactory.create(self.geometry_name, **geometry_args)

        self.interaction_range = kwargs["interaction_range"]
        self.xy_interaction_name = kwargs["xy_interaction_type"] if "xy_interaction_type" in kwargs.keys() else "constant"
        self.xy_interaction_args = InteractionsFactory.extract_args(self.xy_interaction_name, "xy", **kwargs)
        print(self.xy_interaction_args)
        self.xy_interactions = InteractionsFactory.create(
            "xy", 
            self.xy_interaction_name, 
            self.geometry, 
            **self.xy_interaction_args
            )

        self.get_zz_interaction_config(kwargs)

        self.zz_interaction_args = InteractionsFactory.extract_args(self.zz_interaction_name, self.zz_coeffs_type, **kwargs)
        self.zz_interactions = InteractionsFactory.create(
            "zz", 
            self.zz_interaction_name, 
            self.geometry, 
            **self.zz_interaction_args
            )

        self.field_name = kwargs["field_type"] if "field_type" in kwargs.keys() else "constant"
        self.field_args = FieldFactory.extract_args(self.field_name, **kwargs)
        self.field = FieldFactory.create(self.field_name, self.geometry, **self.field_args)

        self.hamiltonian_parameters.construct_coeff_arrays(self.xy_interactions, self.zz_interactions, self.field)
        self.hamiltonian_parameters.validate_coefficient_arrays()

    def get_zz_interaction_config(self, kwarg_dict):
        """
        Sets the interaction configurations for the ZZ term from inputs if available. Defaults to XY case otherwise.

        Args:
            kwarg_dict (dict): contains the interaction_type key value pair

        Returns:
            (str): interaction name
        """

        if "zz_interaction_type" in kwarg_dict.keys():
            self.zz_interaction_name = kwarg_dict["zz_interaction_type"]
            self.zz_coeffs_type = "zz"
        else:
            self.zz_interaction_name = self.xy_interaction_name
            self.zz_coeffs_type = "xy"

    def update_parameters(self, parameter_dict):
        """
        updates a dictionary with system parameters

        Args:
            parameter_dict (dict): parameter set specified as key word arguments

        Returns:
            (dict): updated parameter set
        """

        self.hamiltonian_parameters.update_parameters(parameter_dict)
        self.geometry.update_parameters(parameter_dict)

        self.zz_interactions.update_parameters(parameter_dict)

        return parameter_dict
