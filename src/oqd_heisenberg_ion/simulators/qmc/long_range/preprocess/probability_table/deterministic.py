import os

import numpy as np

from .base import ProbabilityTable


class Deterministic(ProbabilityTable):
    """
    ProbabilityTable subclass for deterministic sampling

    Raises:
        Exception: raised if deterministic probability tables are requested for an unsupported hamiltonian type
        ValueError: raised if the spectrum offset is requested but the hamiltonian type is not consistent with sampling type
    """

    args = {}
    allowed_hamiltonians = {"XY", "fm_heisenberg_afm_Z", "fm_heisenberg_fm_Z"}

    def __init__(self, system):
        """
        validates inputs and constructs the probability tables

        Args:
            system (System): object representing the system to be simulated
        """

        super().__init__(system)

        self.validate_system()

        self.build()

    def validate_system(self):
        """
        validates the system associated with the ProbabilityTable instance

        Raises:
            Exception: if hamiltonian name does not match the allowed hamiltonians for deterministic sampling
        """

        super().validate_system()

        hamiltonian_name = self.system.hamiltonian_parameters.hamiltonian_name

        if hamiltonian_name not in self.allowed_hamiltonians:
            raise Exception(
                "Inconsistent hamiltonian and sampling types. Deterministic probability tables "
                "only support the following types: {}".format(self.allowed_hamiltonians)
            )

    def build(self):
        """
        populates the probability tables for deterministic sampling
        """

        self.compute_max_over_states(self.system.geometry.num_bonds, self.system.xy_interactions.J_ij_vector)
        self.compute_spectrum_offset(self.system.hamiltonian_parameters.hamiltonian_name)

    def compute_spectrum_offset(self, hamiltonian_name):
        """
        computes the spectrum offset associated with the hamiltonian for SSE

        Args:
            hamiltonian_name (str): represents the name of the Hamiltonian

        Raises:
            ValueError: if the spectrum is requested and the Hamiltonian can not be used with deterministic sampling
        """

        if hamiltonian_name == "XY":
            self.spectrum_offset = self.max_diag_norm
        elif hamiltonian_name == "fm_heisenberg_afm_Z" or hamiltonian_name == "fm_heisenberg_fm_Z":
            self.spectrum_offset = 0.5 * self.max_diag_norm
        else:
            raise ValueError(
                "Invalid Hamiltonian type: {} provided for deterministic probability tables. "
                "Allowed types are {}".format(hamiltonian_name, self.allowed_hamiltonians)
            )

    def compute_max_over_states(self, num_bonds, J_ij_vector):
        """
        computes the max norm probability table required for diagonal updates in the two-step method.
        See: https://scipost.org/submissions/2107.00766v1/ and https://journals.aps.org/pre/abstract/10.1103/PhysRevE.66.046701 for details

        Args:
            num_bonds (int): number of interacting bonds in the lattice
            J_ij_vector (numpy.ndarray[float]): num_bonds x 1 array containing the coupling strengths
        """

        max_over_states = np.zeros(num_bonds)
        max_diag_norm = 0.0

        for bond in range(num_bonds):
            J_ij = J_ij_vector[bond]
            max_over_states[bond] = 0.5 * J_ij

            max_diag_norm += 0.5 * J_ij

        self.max_over_states = max_over_states / max_diag_norm
        self.max_diag_norm = max_diag_norm

    def write_to_files(self, out_dir):
        """
        writes the probability tables to files for deterministic sampling

        Args:
            out_dir (str): file path for writing probability tables
        """

        super().write_to_files(out_dir)

        geometry_file_name = os.path.join(self.prob_dir, "geometry.csv")
        max_over_states_file_name = os.path.join(self.prob_dir, "max_over_states.csv")

        geometry_table = self.system.geometry.geometry_table
        num_bonds = self.system.geometry.num_bonds

        np.savetxt(geometry_file_name, geometry_table, delimiter=",", fmt="%d", header="NumBonds={}".format(num_bonds))

        header = "norm={},spectrum_offset={},loop_update_type={}".format(
            self.max_diag_norm, self.spectrum_offset, "deterministic"
        )

        np.savetxt(max_over_states_file_name, self.max_over_states, delimiter=",", header=header)
