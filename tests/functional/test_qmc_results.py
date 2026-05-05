import os

import numpy as np
import pytest

from oqd_heisenberg_ion.common.driver.factory import DriverFactory
from oqd_heisenberg_ion.common.inputs.input_reader import InputReader
from oqd_heisenberg_ion.common.postprocess import utils as post_proc
from oqd_heisenberg_ion.common.preprocessor.factory import PreprocessorFactory


def compute_qmc_energy(preprocessor_qmc):

    qmc_outputs_file = os.path.join(
        preprocessor_qmc.processed_configs[0]["misc"]["run_folder"], "qmc_output/estimators.csv"
    )
    estimator_outputs = np.loadtxt(qmc_outputs_file, delimiter=",", skiprows=2)
    qmc_energy_estimates = estimator_outputs[:, 1]
    energy_stats = post_proc.statistics_binning(qmc_energy_estimates, 2, 100)

    return energy_stats[0], energy_stats[1]


def compute_ed_energy(preprocessor_ed):

    ed_outputs_file = os.path.join(preprocessor_ed.processed_configs[0]["misc"]["run_folder"], "ed_output/energies.csv")
    T = preprocessor_ed.processed_configs[0]["sampling"]["T"]
    ed_outputs = np.loadtxt(ed_outputs_file, delimiter=",", skiprows=1)
    ground_state_energy = post_proc.ed_energy(ed_outputs, T, 1)

    return ground_state_energy


@pytest.mark.parametrize(
    ["hamiltonian_name", "alpha", "h", "loop_type"],
    [
        ["fm_heisenberg_afm_Z", 1.0, 0.0, "deterministic"],
        ["fm_heisenberg_fm_Z", 2.0, 0.0, "deterministic"],
        ["XXZ", 1.5, 0.0, "directed_loop"],
        ["XXZh", 2.5, 1.5, "directed_loop"],
        ["XY", 1.0, 0.0, "deterministic"],
        ["XY", 3.0, 0.0, "deterministic"],
        ["XY", 10.0, 0.0, "deterministic"],
    ],
)
def test_long_range(hamiltonian_name, alpha, h, loop_type, tmp_path):

    input_file = "tests/input_files/long_range.txt"
    inputs = InputReader(input_file_path=input_file)
    inputs.read_inputs_from_file()

    inputs.read_kwarg_inputs(
        hamiltonian_name=hamiltonian_name,
        alpha=alpha,
        h = h,
        loop_type=loop_type,
        equilibration_steps=10000,
        simulation_steps=10000,
        root_folder=tmp_path,
    )

    parameter_set_list = inputs.parameter_set_list
    preprocessor_qmc = PreprocessorFactory.create("long_range_qmc", parameter_set_list)
    driver_inputs = preprocessor_qmc.preprocess()
    driver = DriverFactory.create("long_range_qmc", preprocessor_qmc.simulation_folder, driver_inputs)
    driver.simulate()

    inputs.read_kwarg_inputs(output_folder_name="exact_diagonalization")
    parameter_set_list = inputs.parameter_set_list
    preprocessor_ed = PreprocessorFactory.create("exact_diagonalization", parameter_set_list)
    driver_inputs_ed = preprocessor_ed.preprocess()
    driver_ed = DriverFactory.create("exact_diagonalization", preprocessor_ed.simulation_folder, driver_inputs_ed)
    driver_ed.simulate()

    energy_mean, energy_err = compute_qmc_energy(preprocessor_qmc)
    ed_energy = compute_ed_energy(preprocessor_ed)

    assert energy_mean == pytest.approx(ed_energy, abs=energy_err + 0.1)


@pytest.mark.parametrize(
    ["hamiltonian_name", "J"],
    [
        ["XY", 1.0],
        ["fm_heisenberg_fm_Z", 1.0],
        ["afm_heisenberg_fm_Z", -1.0],
    ],
)
def test_nearest_neighbor(hamiltonian_name, J, tmp_path):

    input_file = "tests/input_files/nearest_neighbor.txt"
    inputs = InputReader(input_file_path=input_file)
    inputs.read_inputs_from_file()

    inputs.read_kwarg_inputs(
        hamiltonian_name=hamiltonian_name, J=J, root_folder=tmp_path, equilibration_steps=10000, simulation_steps=10000
    )
    parameter_set_list = inputs.parameter_set_list
    preprocessor_qmc = PreprocessorFactory.create("nearest_neighbor_qmc", parameter_set_list)
    driver_inputs = preprocessor_qmc.preprocess()
    driver = DriverFactory.create("nearest_neighbor_qmc", preprocessor_qmc.simulation_folder, driver_inputs)
    driver.simulate()

    inputs.read_kwarg_inputs(output_folder_name="exact_diagonalization_nearest_neighbor")
    parameter_set_list = inputs.parameter_set_list
    preprocessor_ed = PreprocessorFactory.create("exact_diagonalization", parameter_set_list)
    driver_inputs_ed = preprocessor_ed.preprocess()
    driver_ed = DriverFactory.create("exact_diagonalization", preprocessor_ed.simulation_folder, driver_inputs_ed)
    driver_ed.simulate()

    energy_mean, energy_err = compute_qmc_energy(preprocessor_qmc)
    ed_energy = compute_ed_energy(preprocessor_ed)

    assert energy_mean == pytest.approx(ed_energy, abs=energy_err + 0.1)


def test_input_matrix(tmp_path):

    input_file = "tests/input_files/long_range.txt"
    inputs = InputReader(input_file_path=input_file)
    inputs.read_inputs_from_file()

    interaction_matrix_file = "tests/input_files/interaction_matrix.csv"

    inputs.read_kwarg_inputs(
        hamiltonian_name="XY",
        interaction_type="matrix_input",
        interaction_matrix_file=interaction_matrix_file,
        loop_type="deterministic",
        equilibration_steps=10000,
        simulation_steps=10000,
        root_folder=tmp_path,
    )

    parameter_set_list = inputs.parameter_set_list
    preprocessor_qmc = PreprocessorFactory.create("long_range_qmc", parameter_set_list)
    driver_inputs = preprocessor_qmc.preprocess()
    driver = DriverFactory.create("long_range_qmc", preprocessor_qmc.simulation_folder, driver_inputs)
    driver.simulate()

    inputs.read_kwarg_inputs(output_folder_name="exact_diagonalization")
    parameter_set_list = inputs.parameter_set_list
    preprocessor_ed = PreprocessorFactory.create("exact_diagonalization", parameter_set_list)
    driver_inputs_ed = preprocessor_ed.preprocess()
    driver_ed = DriverFactory.create("exact_diagonalization", preprocessor_ed.simulation_folder, driver_inputs_ed)
    driver_ed.simulate()

    energy_mean, energy_err = compute_qmc_energy(preprocessor_qmc)
    ed_energy = compute_ed_energy(preprocessor_ed)

    assert energy_mean == pytest.approx(ed_energy, abs=energy_err + 0.1)
