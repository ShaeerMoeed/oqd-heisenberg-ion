import os

import numpy as np

from .base import ProbabilityTable
from .utils import math_utils as mu
from .utils import vertex_utils as vu


class Heatbath(ProbabilityTable):
    """
    ProbabilityTable subclass for heatbath sampling
    See: https://journals.aps.org/pre/abstract/10.1103/PhysRevE.66.046701 for details
    """

    args = {"gamma": float}
    allowed_hamiltonians = {"XXZ", "XXZh", "XY", "fm_heisenberg_fm_Z", "fm_heisenberg_afm_Z", "inhomogenous_XXZh"}

    def __init__(self, system, gamma):
        """
        constructor computes the field contribution per bond, sets member variables and populates the heatbath probability tables

        Args:
            system (System): object representing the system to be simulated
            gamma (float): offset added to weights reduce bounces
        """

        super().__init__(system, gamma=gamma)

        self.gamma = gamma

        self.validate_system()

        self.build()

    def validate_system(self):
        """
        validates the systema associated with the instance of the ProbabilityTable object

        Raises:
            Exception: if, for the specified hamiltonian name, heatbath sampling can not be used
        """

        super().validate_system()

        hamiltonian_name = self.system.hamiltonian_parameters.hamiltonian_name

        if hamiltonian_name not in self.allowed_hamiltonians:
            raise Exception(
                "Inconsistent hamiltonian and sampling types. Heatbath probability tables "
                "only support the following types: {}".format(self.allowed_hamiltonians)
            )

    def build(self):
        """
        populates the heatbath probability tables
        """

        num_bonds = self.system.geometry.num_bonds
        xy_interactions_vector = self.system.hamiltonian_parameters.xy_coeff_array
        zz_interactions_vector = self.system.hamiltonian_parameters.zz_coeff_array
        h_b_vector = self.system.hamiltonian_parameters.z_field_array
        gamma = self.gamma

        self.initialize_tables(num_bonds)
        self.compute_prob_tables_heat_bath(num_bonds, xy_interactions_vector, gamma, h_b_vector, zz_interactions_vector)

        return 0

    def initialize_tables(self, num_bonds):
        """
        initializes the probability tables needed for heatbath sampling

        Args:
            num_bonds (int): number of interacting bonds in the lattice
        """

        self.num_rows = vu.num_vertices * vu.num_legs_indices
        self.heat_bath_prob_table = np.zeros((self.num_rows, num_bonds))

        self.diag_prob_table = np.zeros((vu.num_diagonal_vertices, num_bonds))
        self.max_over_states = np.zeros(num_bonds)
        self.vertex_weights = np.zeros((vu.num_vertices, num_bonds))

        self.spectrum_offset = 0.0
        self.max_diag_norm = 0.0

        return 0

    # Helper function used by compute_prob_tables_heat_bath
    def compute_offset(self, gamma, zz_coeff, h_b):

        if h_b < 0.0:
            raise Exception("h_B needs to be greater than or equal to 0")
        else:
            zz_coeff_over_4 = zz_coeff/4.0

            if zz_coeff_over_4 > h_b:
                offset_b = zz_coeff_over_4

            elif zz_coeff < 0.0:
                offset_b = h_b - zz_coeff_over_4

            else:
                offset_b = h_b

            offset_b += gamma

            return offset_b

    def update_heat_bath_probs(self, bond):

        for vertex_enum in range(vu.num_vertices):
            for l_e in range(vu.num_legs_per_vertex):
                norm = 0.0
                count_invalid_vertices = 0

                for l_x in range(vu.num_legs_per_vertex):
                    composite_leg_index = vu.num_legs_per_vertex * l_e + l_x
                    row_index = vu.num_legs_indices * vertex_enum + composite_leg_index

                    # new_vertex = get_new_vertex(v_map, l_spin, l_e, l_x, vertex_enum)
                    new_vertex = vu.new_vertex_map[vertex_enum, composite_leg_index]

                    if new_vertex == -1:
                        count_invalid_vertices += 1
                        self.heat_bath_prob_table[row_index, bond] = 0.0
                    else:
                        self.heat_bath_prob_table[row_index, bond] = mu.set_probability(
                            self.vertex_weights[new_vertex, bond]
                        )
                        norm += self.vertex_weights[new_vertex, bond]

                self.heat_bath_prob_table[row_index - vu.num_legs_per_vertex + 1 : row_index + 1, bond] /= norm

        return 0

    # Generating this table might be slow for large systems because size grows as N^2.
    # Simpler but slightly slower approach would be to:
    # compute non-zero heat bath probabilities on the fly
    def compute_prob_tables_heat_bath(self, num_bonds, xy_coeff_vector, gamma, h_b_vector, zz_coeff_vector):

        for bond in range(num_bonds):
            xy_coeff = xy_coeff_vector[bond]
            zz_coeff = zz_coeff_vector[bond]
            h_b = h_b_vector[bond]

            vu.set_vertex_weights(self.vertex_weights, bond, zz_coeff, xy_coeff, h_b)

            self.diag_prob_table[:, bond] = self.vertex_weights[0:4, bond]

            offset = self.compute_offset(gamma, zz_coeff, h_b)

            self.vertex_weights[0:4, bond] += offset
            self.spectrum_offset += offset

            self.diag_prob_table[:, bond] += offset

            self.diag_prob_table = mu.enforce_positive(self.diag_prob_table, bond)
            self.vertex_weights = mu.enforce_positive(self.vertex_weights, bond)

            self.max_over_states[bond] = np.max(self.diag_prob_table[:, bond])
            self.diag_prob_table[:, bond] /= self.max_over_states[bond]
            self.max_diag_norm += self.max_over_states[bond]

            self.update_heat_bath_probs(bond)

        self.max_over_states[:] /= self.max_diag_norm

        return 0

    def write_to_files(self, out_dir):
        """
        writes the probability tables to csv files for SSE engine

        Args:
            out_dir (str): directory path for writing probability tables
        """

        super().write_to_files(out_dir)

        geometry_file_name = os.path.join(self.prob_dir, "geometry.csv")
        diag_file_name = os.path.join(self.prob_dir, "diag_probs.csv")
        max_over_states_file_name = os.path.join(self.prob_dir, "max_over_states.csv")
        loop_update_table_file_name = os.path.join(self.prob_dir, "off_diag_table.csv")
        vertex_weights_file_name = os.path.join(self.prob_dir, "vertex_weights.csv")

        geometry_table = self.system.geometry.geometry_table
        num_bonds = self.system.geometry.num_bonds
        np.savetxt(geometry_file_name, geometry_table, delimiter=",", fmt="%d", header="NumBonds={}".format(num_bonds))

        header = "norm={},spectrum_offset={},loop_update_type={}".format(
            self.max_diag_norm, self.spectrum_offset, "heatbath"
        )

        np.savetxt(diag_file_name, self.diag_prob_table, delimiter=",", header=header)
        np.savetxt(vertex_weights_file_name, self.vertex_weights, delimiter=",", header=header)
        np.savetxt(max_over_states_file_name, self.max_over_states, delimiter=",", header=header)
        np.savetxt(loop_update_table_file_name, self.heat_bath_prob_table, delimiter=",", header=header)

        return 0
