from .utils import convert_to_snake_case

input_schema = {
    "hamiltonian_name": {
        "DataType": "categorical",
        "Categories": [
            "fm_heisenberg_afm_Z",
            "XY",
            "fm_heisenberg_fm_Z",
            "XXZ",
            "XXZh",
            "afm_heisenberg_fm_Z",
            "XXZhB",
        ],
        "ParameterType": "system",
        "ConvertValCase": False,
    },
    "loop_type": {
        "DataType": "categorical",
        "Categories": ["deterministic", "heatbath", "directed_loop"],
        "ParameterType": "sampling",
        "ConvertValCase": True,
    },
    "boundary": {
        "DataType": "categorical",
        "Categories": ["periodic", "open"],
        "ParameterType": "system",
        "ConvertValCase": True,
    },
    "interaction_range": {
        "DataType": "categorical",
        "Categories": ["nearest_neighbor", "long_range"],
        "ParameterType": "system",
        "ConvertValCase": True,
    },
    "xy_interaction_type": {
        "DataType": "categorical",
        "Categories": ["power_law", "matrix_input", "vector_input"],
        "ParameterType": "system",
        "ConvertValCase": True,
    },
    "zz_interaction_type": {
        "DataType": "categorical",
        "Categories": ["power_law", "matrix_input", "vector_input"],
        "ParameterType": "system",
        "ConvertValCase": True,
    },
    "spatial_dimension": {
        "DataType": "categorical",
        "Categories": ["1d", "2d"],
        "ParameterType": "system",
        "ConvertValCase": False,
    },
    "lattice_type": {
        "DataType": "categorical",
        "Categories": ["rectangular", "triangular"],
        "ParameterType": "system",
        "ConvertValCase": True,
    },
    "simulator": {
        "DataType": "categorical",
        "Categories": ["long_range_qmc", "nearest_neighbor_qmc", "exact_diagonalization"],
        "ParameterType": "simulation",
        "ConvertValCase": True,
    },
    "N": {"DataType": int, "ParameterType": "system", "ConvertValCase": False},
    "Delta": {"DataType": float, "ParameterType": "system", "ConvertValCase": False},
    "xy_alpha": {"DataType": float, "ParameterType": "system", "ConvertValCase": False},
    "zz_alpha": {"DataType": float, "ParameterType": "system", "ConvertValCase": False},
    "xy_interaction_matrix_file": {"DataType": str, "ParameterType": "system", "ConvertValCase": False},
    "xy_interaction_vector_file": {"DataType": str, "ParameterType": "system", "ConvertValCase": False},
    "zz_interaction_matrix_file": {"DataType": str, "ParameterType": "system", "ConvertValCase": False},
    "zz_interaction_vector_file": {"DataType": str, "ParameterType": "system", "ConvertValCase": False},
    "h": {"DataType": float, "ParameterType": "system", "ConvertValCase": False},
    "z_field_vector_file": {"DataType": str, "ParameterType": "system", "ConvertValCase": False},
    "J": {"DataType": float, "ParameterType": "system", "ConvertValCase": False},
    "T": {"DataType": float, "ParameterType": "sampling", "ConvertValCase": False},
    "B": {"DataType": float, "ParameterType": "system", "ConvertValCase": False},
    "theta": {"DataType": float, "ParameterType": "system", "ConvertValCase": False},
    "equilibration_steps": {"DataType": int, "ParameterType": "simulation", "ConvertValCase": False},
    "simulation_steps": {"DataType": int, "ParameterType": "simulation", "ConvertValCase": False},
    "operator_list_update_multiplier": {"DataType": float, "ParameterType": "simulation", "ConvertValCase": False},
    "gamma": {"DataType": float, "ParameterType": "sampling", "ConvertValCase": False},
    "ksi": {"DataType": float, "ParameterType": "sampling", "ConvertValCase": False},
    "distance_dependent_offset": {"DataType": bool, "ParameterType": "sampling", "ConvertValCase": False},
    "root_folder": {"DataType": str, "ParameterType": "misc", "ConvertValCase": False},
    "bin_folder": {"DataType": str, "ParameterType": "misc", "ConvertValCase": False},
    "cpp_source_folder": {"DataType": str, "ParameterType": "misc", "ConvertValCase": False},
    "julia_path": {"DataType": str, "ParameterType": "misc", "ConvertValCase": False},
    "uuid": {"DataType": str, "ParameterType": "misc", "ConvertValCase": False},
    "output_folder_name": {"DataType": str, "ParameterType": "misc", "ConvertValCase": False},
    "number_of_threads": {"DataType": int, "ParameterType": "misc", "ConvertValCase": False},
    "track_spin_configurations": {"DataType": bool, "ParameterType": "simulation", "ConvertValCase": False},
    "write_final_spin_configuration": {"DataType": bool, "ParameterType": "simulation", "ConvertValCase": False},
    "initial_configuration_index": {"DataType": int, "ParameterType": "simulation", "ConvertValCase": False},
    "initial_configuration_file_path": {"DataType": str, "ParameterType": "simulation", "ConvertValCase": False},
    "initial_operator_list_size": {"DataType": int, "ParameterType": "simulation", "ConvertValCase": False},
    "initial_config_seed": {"DataType": int, "ParameterType": "simulation", "ConvertValCase": False},
    "disconnected_spin_flip_seed": {"DataType": int, "ParameterType": "simulation", "ConvertValCase": False},
    "off_diagonal_update_seed": {"DataType": int, "ParameterType": "simulation", "ConvertValCase": False},
    "metropolis_insert_seed": {"DataType": int, "ParameterType": "simulation", "ConvertValCase": False},
    "metropolis_remove_seed": {"DataType": int, "ParameterType": "simulation", "ConvertValCase": False},
    "metropolis_bond_generator_seed": {"DataType": int, "ParameterType": "simulation", "ConvertValCase": False},
    "diagonal_update_seed": {"DataType": int, "ParameterType": "simulation", "ConvertValCase": False},
    "exit_leg_seed": {"DataType": int, "ParameterType": "simulation", "ConvertValCase": False},
}


class InputParser:
    """
    Parses the provided list of inputs for a single parameter set.
    Converts strings to required types using the data schema and builds the input configurations for the driver
    """

    def __init__(self, **config_settings):
        """
        constructs an instance of InputParser. Builds the input_config settings for the driver layer

        Args:
            **config_settings (dict): key word arguments specifying a single parameter set with str values
        """

        self.dtype_parsers = {
            str: self.extract_string,
            int: self.extract_integer,
            float: self.extract_float,
            bool: self.extract_bool,
            "categorical": self.extract_categorical,
        }

        self.input_schema = input_schema

        self.simulation_config = {}

        self.build_input_config(config_settings)

    def extract_integer(self, config_settings, key):
        """
        converts str to int for parameter specified by key

        Args:
            config_settings (dict): contains the key value pairs that define the simulation parameter set
            key (str): parameter key

        Returns:
            (int): parameter value
        """

        return int(config_settings[key])

    def extract_float(self, config_settings, key):
        """
        converts str to float for parameter specified by key

        Args:
            config_settings (dict): contains the key value pairs that define the simulation parameter set
            key (str): parameter key

        Returns:
            (float): parameter value
        """

        return float(config_settings[key])

    def extract_string(self, config_settings, key):
        """
        extracts str parameter specified by key

        Args:
            config_settings (dict): contains the key value pairs that define the simulation parameter set
            key (str): parameter key

        Returns:
            (str): parameter value
        """

        return config_settings[key]

    def extract_bool(self, config_settings, key):
        """
        converts str to bool for parameter specified by key

        Args:
            config_settings (dict): contains the key value pairs that define the simulation parameter set
            key (str): parameter key

        Returns:
            (bool): parameter value
        """

        if isinstance(config_settings[key], bool):
            return config_settings[key]
        else:
            if config_settings[key].capitalize() == "True":
                return True
            elif config_settings[key].capitalize() == "False":
                return False
            else:
                raise ValueError(
                    f"Unrecognized entry for key {key}, value provided: {config_settings[key]}. "
                    "Allowed values are 'True' or 'False'\n"
                )

    def extract_categorical(self, config_settings, key):
        """
        extracts categorical parameter specified by key. Checks if the provided value is allowed by schema

        Args:
            config_settings (dict): contains the key value pairs that define the simulation parameter set
            key (str): parameter key

        Raises:
            ValueError: if value associated with key is unrecognized by schema

        Returns:
            (str): parameter value
        """

        val = config_settings[key]

        allowed_vals = self.input_schema[key]["Categories"]

        if val in allowed_vals:
            return convert_to_snake_case(val) if self.input_schema[key]["ConvertValCase"] else val
        else:
            raise ValueError(
                f"Unrecognized input for key: {key}. Value provided: {val}. Allowed values: {allowed_vals}"
            )

    def build_input_config(self, config_settings):
        """
        builds the entire input configuration corresponding to the parameter set with correct types

        Args:
            config_settings (dict): key word arguments specifying a single parameter set with str values

        """

        for key in config_settings.keys():
            dtype = self.input_schema[key]["DataType"]
            param_type = self.input_schema[key]["ParameterType"]

            parser = self.dtype_parsers[dtype]
            val = parser(config_settings, key)

            if param_type not in self.simulation_config:
                self.simulation_config[param_type] = {}

            self.simulation_config[param_type][key] = val
