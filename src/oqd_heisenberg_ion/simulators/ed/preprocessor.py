import os
import shutil as sh

from oqd_heisenberg_ion.common.inputs.input_parser import InputParser
from oqd_heisenberg_ion.common.preprocessor.base import Preprocessor

from ..preprocess.system.base import System


class ExactDiagonalization(Preprocessor):
    """
    Preprocessor subclass for Exact Diagonalization. Validates inputs, configures the parameter sets and writes the engine input file
    """

    def __init__(self, parameter_set_list):
        """
        constructor for the ED Preprocessor.

        Args:
            parameter_set_list (list[dict]): list of parameter sets specified with parameter set defined by a dict
        """

        super().__init__(parameter_set_list)

        # self.driver_inputs = []

    def preprocess(self):
        """
        Preprocesses the simulation parameter sets, and validates root output folder and uuids.
        Also populates the ED Driver inputs.
        """

        self.check_single_input("root_folder")

        self.simulation_folder = self.create_output_folder()

        self.extract_cli_requirements()

        self.check_unique_uuids()

        self.configure_simulation()

        return self.driver_inputs

    def configure_simulation(self):
        """
        Configures all the parameter sets sequentially. Then writes the input file for the engine.
        """

        for i in range(self.num_parameter_sets):
            parameter_set = self.parameter_set_list[i]

            self.configure_parameter_set(parameter_set)

        self.write_input_file()

    def configure_parameter_set(self, parameter_args):
        """
        Configures a single parameter set.
        Determines the run_id, creates the parameter set output folder and parses the str inputs.
        Also defines the system and appends any required parameters to the parameter set.

        Args:
            parameter_args (dict): Parameter set arguments specified as key value pairs. Single value for each key specified as a str
        """

        input_config = InputParser(**parameter_args)
        system_args = input_config.simulation_config["system"]

        misc_args = input_config.simulation_config["misc"]
        run_id = self.get_run_id(misc_args)
        misc_args["uuid"] = run_id
        misc_args["output_folder_name"] = self.output_folder_name

        misc_args["simulation_folder"] = self.simulation_folder

        run_folder = self.create_run_folder(misc_args)
        misc_args["run_folder"] = run_folder

        system = System(**system_args)
        system_args = system.update_parameters(system_args)

        if system.interaction_range == "long_range":
            J_ij_file_path = os.path.join(run_folder, "J_ij_file.csv")
            system.xy_interactions.write_to_file(J_ij_file_path)
            misc_args["J_ij_file"] = J_ij_file_path

        self.processed_configs.append(input_config.simulation_config)

    def extract_cli_requirements(self):
        """
        Extracts the Julia path if specified as a user input. If not provided, attempts to find the Julia path and throws an error if unsucessful.
        If the path is available, populates the ED Driver inputs
        """

        julia_path = self.extract_optional_input("julia_path", True)
        if julia_path is None:
            julia_path = sh.which("julia")
            if julia_path is None:
                raise Exception("No Julia path provided\n")

        self.driver_inputs = {"julia_path": julia_path}
