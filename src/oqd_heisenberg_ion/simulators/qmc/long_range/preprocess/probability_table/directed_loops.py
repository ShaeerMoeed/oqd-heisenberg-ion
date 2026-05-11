import os

import numpy as np

from .base import ProbabilityTable
from .utils import math_utils as mu
from .utils import vertex_utils as vu


class DirectedLoops(ProbabilityTable):
    """
    ProbabilityTable subclass for directed loops sampling
    See: https://journals.aps.org/pre/abstract/10.1103/PhysRevE.66.046701 for details
    """

    args = {"gamma": float, "ksi": float, "distance_dependent_offset": bool}
    allowed_hamiltonians = {"XXZ", "XXZh", "XY", "fm_heisenberg_fm_Z", "fm_heisenberg_afm_Z", "inhomogenous_XXZh"}

    def __init__(self, system, gamma, ksi, distance_dependent_offset):
        """
        constructor computes the field contribution per bond, sets member variables and populates the directed loop probability tables

        Args:
            system (System): object representing the system to be simulated
            gamma (float): offset added to weights reduce bounces
            ksi (float): offset added in the presence of a large field to reduce bounces
            distance_dependent_offset (bool): determines whether a distance dependant offset should be used
        """

        super().__init__(system, gamma=gamma, ksi=ksi, distance_dependent_offset=distance_dependent_offset)

        self.gamma = gamma
        self.ksi = ksi
        self.distance_dependent_offset = distance_dependent_offset

        print("Test")

        self.build()

    def validate_system(self):
        """
        validates the systema associated with the instance of the ProbabilityTable object

        Raises:
            Exception: if, for the specified hamiltonian name, directed loop sampling can not be used
        """

        super().validate_system()

        hamiltonian_name = self.system.hamiltonian_parameter.hamiltonian_name

        if hamiltonian_name not in self.allowed_hamiltonians:
            raise Exception(
                "Inconsistent hamiltonian and sampling types. Directed loop probability tables "
                "only support the following types: {}".format(self.allowed_hamiltonians)
            )

    def build(self):
        """
        populates the directed loop probability tables
        """

        num_bonds = self.system.geometry.num_bonds

        self.initialize_tables(num_bonds)
        self.set_vertex_enum_transition_weights_map()

        xy_coeff_vector = self.system.hamiltonian_parameters.xy_coeff_array
        zz_coeff_vector = self.system.hamiltonian_parameters.zz_coeff_array
        h_b_vector = self.system.hamiltonian_parameters.z_field_array

        gamma = self.gamma
        ksi = self.ksi
        distance_dependent_offset = self.distance_dependent_offset

        self.compute_prob_tables_directed_loops(
            num_bonds, xy_coeff_vector, gamma, h_b_vector, zz_coeff_vector, ksi, distance_dependent_offset
        )

    def initialize_tables(self, num_bonds):
        """
        initializes the probability tables needed for directed loops sampling

        Args:
            num_bonds (int): number of interacting bonds in the lattice
        """

        num_rows = vu.num_vertices * vu.num_legs_indices

        self.directed_loop_prob_table = np.zeros((num_rows, num_bonds))
        self.diag_prob_table = np.zeros((vu.num_diagonal_vertices, num_bonds))
        self.max_over_states = np.zeros(num_bonds)
        self.vertex_weights = np.zeros((vu.num_vertices, num_bonds))

        self.spectrum_offset = 0.0
        self.max_diag_norm = 0.0

    def set_vertex_enum_transition_weights_map(self):

        self.vertex_weight_label_map = {}

        self.num_vertex_enums = 4
        self.directed_loop_vertex_enums = ["0", "1", "5", "3"]
        self.vertex_enum_weight_list_counts = [1, 2, 2, 1]

        vertex_enum = "0"
        exit_leg_weights_le_0_v_0 = ["b_3", None, "c", "b"]  # l_e = 0
        self.vertex_weight_label_map[vertex_enum] = [exit_leg_weights_le_0_v_0]

        vertex_enum = "1"
        exit_leg_weights_le_0_v_1 = ["b_2_p", "a_p", "c_p", None]  # l_e = 0
        exit_leg_weights_le_1_v_1 = ["a", "b_2", None, "c"]  # l_e = 1
        self.vertex_weight_label_map[vertex_enum] = [exit_leg_weights_le_0_v_1, exit_leg_weights_le_1_v_1]

        vertex_enum = "5"
        exit_leg_weights_le_0_v_5 = ["b_1", "a", None, "b"]  # l_e = 0
        exit_leg_weights_le_1_v_5 = ["a_p", "b_1_p", "b_p", None]  # l_e = 1
        self.vertex_weight_label_map[vertex_enum] = [exit_leg_weights_le_0_v_5, exit_leg_weights_le_1_v_5]

        vertex_enum = "3"
        exit_leg_weights_le_0_v_3 = ["b_3_p", None, "c_p", "b_p"]  # l_e = 0
        self.vertex_weight_label_map[vertex_enum] = [exit_leg_weights_le_0_v_3]

    def update_directed_loop_probs(self, vertex_enum, l_e, l_x, bond, transition_weight):

        init_composite_leg_index, init_row_index = self.get_composite_row_prob_index(vertex_enum, l_e, l_x)
        normalization = self.vertex_weights[vertex_enum, bond]
        if normalization != 0.0:
            self.directed_loop_prob_table[init_row_index, bond] = mu.set_probability(transition_weight / normalization)

        new_vertex_enum, new_l_e, new_l_x = self.get_symmetric_indices(vertex_enum, l_e, l_x, vu.vertical_swap_mapping)
        new_composite_leg_index, new_row_index = self.get_composite_row_prob_index(new_vertex_enum, new_l_e, new_l_x)
        self.directed_loop_prob_table[new_row_index, bond] = self.directed_loop_prob_table[init_row_index, bond]

        new_vertex_enum, new_l_e, new_l_x = self.get_symmetric_indices(
            vertex_enum, l_e, l_x, vu.horizontal_swap_mapping
        )
        new_composite_leg_index, new_row_index = self.get_composite_row_prob_index(new_vertex_enum, new_l_e, new_l_x)
        self.directed_loop_prob_table[new_row_index, bond] = self.directed_loop_prob_table[init_row_index, bond]

        new_vertex_enum, new_l_e, new_l_x = self.get_symmetric_indices(vertex_enum, l_e, l_x, vu.composed_swaps_mapping)
        new_composite_leg_index, new_row_index = self.get_composite_row_prob_index(new_vertex_enum, new_l_e, new_l_x)
        self.directed_loop_prob_table[new_row_index, bond] = self.directed_loop_prob_table[init_row_index, bond]

    def get_composite_row_prob_index(self, vertex_enum, entrance_leg_enum, exit_leg_enum):

        composite_leg_index = vu.num_legs_per_vertex * entrance_leg_enum + exit_leg_enum
        row_index = vu.num_legs_indices * vertex_enum + composite_leg_index

        return composite_leg_index, row_index

    def get_symmetric_indices(self, vertex_enum, entrance_leg_enum, exit_leg_enum, symmetry_leg_mapping):

        init_spin_tuple = vu.leg_spin[vertex_enum]
        new_spin_tuple = (
            init_spin_tuple[symmetry_leg_mapping[0]],
            init_spin_tuple[symmetry_leg_mapping[1]],
            init_spin_tuple[symmetry_leg_mapping[2]],
            init_spin_tuple[symmetry_leg_mapping[3]],
        )
        new_vertex_enum = vu.vertex_map[new_spin_tuple]

        new_entrance_leg_enum = symmetry_leg_mapping[entrance_leg_enum]
        new_exit_leg_enum = symmetry_leg_mapping[exit_leg_enum]

        return new_vertex_enum, new_entrance_leg_enum, new_exit_leg_enum

    def update_directed_loop_table(self, bond, transition_weights):

        for i in range(self.num_vertex_enums):
            vertex_enum = self.directed_loop_vertex_enums[i]
            num_unique_entrance_legs = self.vertex_enum_weight_list_counts[i]

            for l_e in range(num_unique_entrance_legs):
                exit_leg_weight_labels = self.vertex_weight_label_map[vertex_enum][l_e]

                for l_x in range(vu.num_legs_per_vertex):
                    l_x_weight_label = exit_leg_weight_labels[l_x]
                    l_x_weight = transition_weights[l_x_weight_label] if l_x_weight_label is not None else 0.0

                    self.update_directed_loop_probs(int(vertex_enum), l_e, l_x, bond, l_x_weight)

    def compute_prob_tables_directed_loops(
        self, num_bonds, xy_coeff_vector, gamma, h_b_vector, zz_coeff_vector, ksi, distance_dependent_offset
    ):

        self.transition_weights_calculator = LoopTransitionWeights(gamma, ksi, distance_dependent_offset)

        for bond in range(num_bonds):
            xy_coeff = xy_coeff_vector[bond]
            zz_coeff = zz_coeff_vector[bond]
            h_b = h_b_vector[bond]

            vu.set_vertex_weights(self.vertex_weights, bond, zz_coeff, xy_coeff, h_b)

            self.diag_prob_table[:, bond] = self.vertex_weights[0:4, bond]

            self.transition_weights_calculator.compute_transition_weights(xy_coeff, zz_coeff, h_b)

            transition_weights = self.transition_weights_calculator.transition_weight_container
            offset = self.transition_weights_calculator.offset_b

            self.vertex_weights[0:4, bond] += offset
            self.spectrum_offset += offset

            self.diag_prob_table[:, bond] += offset

            self.diag_prob_table = mu.enforce_positive(self.diag_prob_table, bond)
            self.vertex_weights = mu.enforce_positive(self.vertex_weights, bond)

            self.update_directed_loop_table(bond, transition_weights)

            self.max_over_states[bond] = np.max(self.diag_prob_table[:, bond])
            self.diag_prob_table[:, bond] /= self.max_over_states[bond]
            self.max_diag_norm += self.max_over_states[bond]

        self.max_over_states[:] /= self.max_diag_norm

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
            self.max_diag_norm, self.spectrum_offset, "DirectedLoops"
        )

        np.savetxt(diag_file_name, self.diag_prob_table, delimiter=",", header=header)
        np.savetxt(vertex_weights_file_name, self.vertex_weights, delimiter=",", header=header)
        np.savetxt(max_over_states_file_name, self.max_over_states, delimiter=",", header=header)
        np.savetxt(loop_update_table_file_name, self.directed_loop_prob_table, delimiter=",", header=header)


class LoopTransitionWeights:
    """
    contains the logic for computing the transition weights in the directed loop SSE method
    See https://journals.aps.org/pre/abstract/10.1103/PhysRevE.66.046701 for details
    """

    def __init__(self, gamma, ksi, distance_dependent_offset):

        keys = ["a", "b", "c", "a_p", "b_p", "c_p", "b_1", "b_2", "b_3", "b_1_p", "b_2_p", "b_3_p"]
        self.transition_weight_container = {key: None for key in keys}

        self.gamma = gamma
        self.ksi = ksi
        self.distance_dependent_offset = distance_dependent_offset

    def populate_unprimed_transition_weights(self, a, b, c):

        self.transition_weight_container["a"] = a
        self.transition_weight_container["b"] = b
        self.transition_weight_container["c"] = c

    def populate_primed_transition_weights(self, a_p, b_p, c_p):

        self.transition_weight_container["a_p"] = a_p
        self.transition_weight_container["b_p"] = b_p
        self.transition_weight_container["c_p"] = c_p

    def populate_unprimed_bounce_weights(self, b_1, b_2, b_3):

        self.transition_weight_container["b_1"] = b_1
        self.transition_weight_container["b_2"] = b_2
        self.transition_weight_container["b_3"] = b_3

    def populate_primed_bounce_weights(self, b_1_p, b_2_p, b_3_p):

        self.transition_weight_container["b_1_p"] = b_1_p
        self.transition_weight_container["b_2_p"] = b_2_p
        self.transition_weight_container["b_3_p"] = b_3_p

    def swap_primed_unprimed_weights(self):

        a = self.transition_weight_container["a"]
        b = self.transition_weight_container["b"]
        c = self.transition_weight_container["c"]

        a_p = self.transition_weight_container["a_p"]
        b_p = self.transition_weight_container["b_p"]
        c_p = self.transition_weight_container["c_p"]

        self.populate_unprimed_transition_weights(a_p, b_p, c_p)
        self.populate_primed_transition_weights(a, b, c)

        b_1 = self.transition_weight_container["b_1"]
        b_2 = self.transition_weight_container["b_2"]
        b_3 = self.transition_weight_container["b_3"]

        b_1_p = self.transition_weight_container["b_1_p"]
        b_2_p = self.transition_weight_container["b_2_p"]
        b_3_p = self.transition_weight_container["b_3_p"]

        self.populate_unprimed_bounce_weights(b_1_p, b_2_p, b_3_p)
        self.populate_primed_bounce_weights(b_1, b_2, b_3)

    def tranisiton_weights_small_field(self, zz_coeff_over_four, zz_p_xy_over_2, zz_m_xy_over_2, h_b):

        self.offset_b = zz_coeff_over_four

        if h_b >= zz_m_xy_over_2:
            b_3_p = 0.0
            epsilon = -zz_m_xy_over_2 / 2.0 + h_b / 2.0 + self.gamma
        else:
            b_3_p = zz_m_xy_over_2 - h_b + self.ksi
            epsilon = self.gamma

        if h_b <= -zz_m_xy_over_2:
            b_3 = 0.0
        else:
            b_3 = zz_m_xy_over_2 + h_b + self.ksi

        a_p = -zz_m_xy_over_2 / 2.0 + h_b / 2.0 + b_3_p / 2.0
        b_p = zz_p_xy_over_2 / 2.0 - h_b / 2.0 - b_3_p / 2.0
        c_p = zz_m_xy_over_2 / 2.0 + epsilon - h_b / 2.0 - b_3_p / 2.0

        a = -zz_m_xy_over_2 / 2.0 - h_b / 2.0 + b_3 / 2.0
        b = zz_p_xy_over_2 / 2.0 + h_b / 2.0 - b_3 / 2.0
        c = epsilon + zz_m_xy_over_2 / 2.0 + h_b / 2.0 - b_3 / 2.0

        b_1_p = 0.0
        b_2_p = 0.0
        b_1 = 0.0
        b_2 = 0.0

        self.offset_b += epsilon
        self.epsilon = epsilon

        self.populate_unprimed_transition_weights(a, b, c)
        self.populate_primed_transition_weights(a_p, b_p, c_p)
        self.populate_unprimed_bounce_weights(b_1, b_2, b_3)
        self.populate_primed_bounce_weights(b_1_p, b_2_p, b_3_p)

    def transition_weights_negative_Delta(self, zz_coeff_over_four, zz_p_xy_over_2, zz_m_xy_over_2, xy_coeff, h_b):

        if h_b == 0.0:
            self.offset_b = -zz_coeff_over_four
            if self.Delta <= -1.0:
                if self.distance_dependent_offset:
                    epsilon = self.gamma - zz_coeff_over_four/2.5
                    c_p = self.gamma - zz_coeff_over_four/2.5
                else:
                    epsilon = self.gamma
                    c_p = self.gamma
                c = c_p
                a_p = (1.0 / 2.0) * xy_coeff
                b_p = 0.0
                a = a_p
                b = b_p
                b_2_p = -zz_p_xy_over_2
                b_2 = b_2_p
                b_3_p = 0.0
                b_1_p = 0.0
                b_1 = b_1_p
                b_3 = b_3_p
            else:
                b_2_p = 0.0
                b_p = zz_p_xy_over_2
                a_p = -zz_m_xy_over_2
                c_p = self.gamma
                epsilon = zz_p_xy_over_2/2.0 + self.gamma
                c = c_p
                a = a_p
                b = b_p
                b_1_p = 0.0
                b_3_p = 0.0
                b_2 = b_2_p
                b_1 = b_1_p
                b_3 = b_3_p
        else:
            self.offset_b = h_b - zz_coeff_over_four

            if h_b <= zz_p_xy_over_2:
                b_2_p = 0.0
                epsilon = zz_p_xy_over_2 / 2.0 - h_b / 2.0 + self.gamma
            else:
                b_2_p = h_b - zz_p_xy_over_2 + self.ksi
                epsilon = self.gamma

            if h_b <= -zz_p_xy_over_2:
                b_2 = -h_b - zz_p_xy_over_2 + self.ksi
            else:
                b_2 = 0.0

            if h_b <= -zz_m_xy_over_2:
                b_3 = 0.0
            else:
                b_3 = h_b + zz_m_xy_over_2 + self.ksi

            a_p = -zz_m_xy_over_2 / 2.0 + h_b / 2.0 - b_2_p / 2.0
            b_p = zz_p_xy_over_2 / 2.0 - h_b / 2.0 + b_2_p / 2.0
            c_p = epsilon - zz_p_xy_over_2 / 2.0 + h_b / 2.0 - b_2_p / 2.0

            a = -zz_m_xy_over_2 / 2.0 - h_b / 2.0 + b_3 / 2.0 - b_2 / 2.0
            b = zz_p_xy_over_2 / 2.0 + h_b / 2.0 + b_2 / 2.0 - b_3 / 2.0
            c = 3.0 * h_b / 2.0 + epsilon - zz_p_xy_over_2 / 2.0 - b_2 / 2.0 - b_3 / 2.0

            b_1 = 0.0
            b_1_p = 0.0
            b_3_p = 0.0

        self.offset_b += epsilon
        self.epsilon = epsilon

        self.populate_unprimed_transition_weights(a, b, c)
        self.populate_primed_transition_weights(a_p, b_p, c_p)
        self.populate_unprimed_bounce_weights(b_1, b_2, b_3)
        self.populate_primed_bounce_weights(b_1_p, b_2_p, b_3_p)

    def transition_weights_large_field(self, zz_p_xy_over_2, zz_m_xy_over_2, xy_coeff, h_b):

        self.offset_b = h_b
        xy_coeff_over_four = xy_coeff/4.0

        if h_b <= zz_m_xy_over_2:
            b_3_p = zz_m_xy_over_2 - h_b + self.ksi
            epsilon = self.gamma
        else:
            b_3_p = 0.0
            if h_b <= zz_p_xy_over_2 and h_b <= 2.0 * xy_coeff_over_four:
                epsilon = xy_coeff_over_four - h_b / 2.0 + self.gamma
            else:
                epsilon = self.gamma

        if h_b <= zz_p_xy_over_2:
            b_2_p = 0.0
        else:
            b_2_p = h_b - zz_p_xy_over_2 + self.ksi

        if h_b < -zz_m_xy_over_2:
            b_3 = 0
        else:
            b_3 = h_b + zz_m_xy_over_2 + self.ksi

        a_p = -zz_m_xy_over_2 / 2.0 + h_b / 2.0 + b_3_p / 2.0 - b_2_p / 2.0
        b_p = zz_p_xy_over_2 / 2.0 - h_b / 2.0 - b_3_p / 2.0 + b_2_p / 2.0
        c_p = epsilon - xy_coeff_over_four + h_b / 2.0 - b_3_p / 2.0 - b_2_p / 2.0

        a = -zz_m_xy_over_2 / 2.0 - h_b / 2.0 + b_3 / 2.0
        b = zz_p_xy_over_2 / 2.0 + h_b / 2.0 - b_3 / 2.0
        c = epsilon - xy_coeff_over_four + 3.0 * h_b / 2.0 - b_3 / 2.0

        b_1 = 0.0
        b_2 = 0.0
        b_1_p = 0.0

        self.offset_b += epsilon
        self.epsilon = epsilon

        self.populate_unprimed_transition_weights(a, b, c)
        self.populate_primed_transition_weights(a_p, b_p, c_p)
        self.populate_unprimed_bounce_weights(b_1, b_2, b_3)
        self.populate_primed_bounce_weights(b_1_p, b_2_p, b_3_p)

    def compute_transition_weights(self, xy_coeff, zz_coeff, h_b):

        zz_coeff_over_four = zz_coeff/4.0

        zz_p_xy_over_2 = zz_coeff/2.0 + xy_coeff/2.0
        zz_m_xy_over_2 = zz_coeff/2.0 - xy_coeff/2.0

        if zz_coeff_over_four > np.abs(h_b):
            self.tranisiton_weights_small_field(zz_coeff_over_four, zz_p_xy_over_2, zz_m_xy_over_2, np.abs(h_b))
        elif zz_coeff < 0.0:
            self.transition_weights_negative_Delta(zz_coeff_over_four, zz_p_xy_over_2, zz_m_xy_over_2, xy_coeff, np.abs(h_b))
        else:
            self.transition_weights_large_field(zz_p_xy_over_2, zz_m_xy_over_2, xy_coeff, np.abs(h_b))

        if h_b < 0.0:
            self.swap_primed_unprimed_weights()
