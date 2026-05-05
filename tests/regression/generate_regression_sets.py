import json
import os

import numpy as np

from oqd_heisenberg_ion.common.driver.factory import DriverFactory
from oqd_heisenberg_ion.common.inputs.input_reader import InputReader
from oqd_heisenberg_ion.common.preprocessor.factory import PreprocessorFactory

data_for_json = []


def compute_qmc_stats(qmc_outputs_file):
    estimator_outputs = np.loadtxt(qmc_outputs_file, delimiter=",", skiprows=2)
    return np.mean(estimator_outputs, axis=0), np.std(estimator_outputs, axis=0)


# Long range QMC
#names = ["fm_heisenberg_afm_Z", "fm_heisenberg_fm_Z", "XXZ", "XXZh", "XXZh", "XY", "XY", "XY"]
#alphas = [1.0, 2.0, 1.5, 2.5, 3.0, 1.0, 3.0, 10.0]
#h_list = [0.0, 0.0, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0]
#loop_types = ["deterministic", "deterministic", "directed_loop", "directed_loop", "heatbath", "deterministic", "deterministic", "deterministic"]

names = ["XXZ", "XXZh"]
alphas = [2.5, 2.5]
h_list = [0.0, 0.0]
loop_types = ["directed_loop", "directed_loop"]

uuids = [names[i] + "_alpha_" + str(alphas[i]) for i in range(len(names))]
input_file = "tests/input_files/long_range.txt"
inputs = InputReader(input_file_path=input_file)
inputs.read_inputs_from_file()
inputs.read_kwarg_inputs(hamiltonian_name=names, alpha=alphas, uuid=uuids, loop_type=loop_types, h=h_list)
parameter_set_list = inputs.parameter_set_list
preprocessor = PreprocessorFactory.create("long_range_qmc", parameter_set_list)
driver_inputs = preprocessor.preprocess()
driver = DriverFactory.create("long_range_qmc", preprocessor.simulation_folder, driver_inputs)
driver.simulate()

for i in range(len(names)):
    qmc_outputs_file = os.path.join(
        preprocessor.processed_configs[i]["misc"]["run_folder"], "qmc_output/estimators.csv"
    )
    estimator_outputs = compute_qmc_stats(qmc_outputs_file)
    stats_dict = {}
    stats_dict["name"] = names[i]
    stats_dict["alpha"] = alphas[i]
    stats_dict["energy_mean"] = estimator_outputs[0][1]
    stats_dict["energy_std"] = estimator_outputs[1][1]
    stats_dict["magnetization_mean"] = estimator_outputs[0][2]
    stats_dict["magnetization_std"] = estimator_outputs[1][2]
    stats_dict["stiffness_mean"] = estimator_outputs[0][3]
    stats_dict["stiffness_std"] = estimator_outputs[1][3]

    data_for_json.append(stats_dict)

# Nearest Neighbor QMC
names = ["XY", "fm_heisenberg_fm_Z", "afm_heisenberg_fm_Z"]
uuids = names
J_list = [1.0, 1.0, -1.0]
input_file = "tests/input_files/nearest_neighbor.txt"
inputs = InputReader(input_file_path=input_file)
inputs.read_inputs_from_file()
inputs.read_kwarg_inputs(hamiltonian_name=names, uuid=uuids, J=J_list)
parameter_set_list = inputs.parameter_set_list
preprocessor = PreprocessorFactory.create("nearest_neighbor_qmc", parameter_set_list)
driver_inputs = preprocessor.preprocess()
driver = DriverFactory.create("nearest_neighbor_qmc", preprocessor.simulation_folder, driver_inputs)
driver.simulate()

for i in range(len(names)):
    qmc_outputs_file = os.path.join(
        preprocessor.processed_configs[i]["misc"]["run_folder"], "qmc_output/estimators.csv"
    )
    estimator_outputs = compute_qmc_stats(qmc_outputs_file)
    stats_dict = {}
    stats_dict["name"] = names[i]
    stats_dict["J"] = J_list[i]
    stats_dict["energy_mean"] = estimator_outputs[0][1]
    stats_dict["energy_std"] = estimator_outputs[1][1]
    stats_dict["magnetization_mean"] = estimator_outputs[0][2]
    stats_dict["magnetization_std"] = estimator_outputs[1][2]

    data_for_json.append(stats_dict)

json_file_path = os.path.abspath("tests/regression/regression_results.json")
with open(json_file_path, "w") as json_file:
    json.dump(data_for_json, json_file, indent=2)
