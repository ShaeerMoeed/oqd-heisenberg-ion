import numpy as np

num_vertices = 6
num_legs_per_vertex = 4
num_diagonal_vertices = 4
num_legs_indices = num_legs_per_vertex**2

# Mapping from vertex to integer.
# The 4 integers that define a vertex specify the spin at each leg (up or down)
vertex_map = {
    (1, 1, 1, 1): 0,
    (1, -1, 1, -1): 1,
    (-1, 1, -1, 1): 2,
    (-1, -1, -1, -1): 3,
    (1, -1, -1, 1): 4,
    (-1, 1, 1, -1): 5,
}

operator_type = (0, 0, 0, 0, 1, 1)  # 0 -> diagonal operator, 1 -> off-diagonal. Index defines vertex type

# Inverse map of vertex_map
leg_spin = [(1, 1, 1, 1), (1, -1, 1, -1), (-1, 1, -1, 1), (-1, -1, -1, -1), (1, -1, -1, 1), (-1, 1, 1, -1)]
# <1|S_z|1> = 1/2, <-1|S_z|-1> = -1/2

# mapping[x] -> y => vertex defined by x maps to y under swap vertex operation
vertical_swap_mapping = [2, 3, 0, 1]
horizontal_swap_mapping = [1, 0, 3, 2]
composed_swaps_mapping = [3, 2, 1, 0]


def get_new_vertex(l_e, l_x, vertex_type):

    if l_e == l_x:
        return vertex_type
    else:
        spin_configs = list(leg_spin[vertex_type])
        spin_configs[l_e] = -spin_configs[l_e]
        spin_configs[l_x] = -spin_configs[l_x]

        spin_configs_tuple = tuple(spin_configs)

        if spin_configs_tuple in leg_spin:
            return vertex_map[tuple(spin_configs)]
        else:
            return -1


def generate_new_vertex_type_array():

    num_cols = num_legs_per_vertex**2
    new_vertex_types = np.zeros((num_vertices, num_cols), dtype=int)
    for i in range(num_vertices):
        for j in range(num_cols):
            l_e = j // num_legs_per_vertex
            l_x = j % num_legs_per_vertex
            new_vertex_types[i, j] = get_new_vertex(l_e, l_x, i)

    return new_vertex_types


new_vertex_map = generate_new_vertex_type_array()  # all possible transitions for all vertices


def set_vertex_weights(vertex_weights, bond, zz_coeff, xy_coeff, h_b):

    vertex_weights[0, bond] = zz_coeff/4.0 + h_b
    vertex_weights[1, bond] = -zz_coeff/4.0
    vertex_weights[2, bond] = vertex_weights[1, bond]
    vertex_weights[3, bond] = zz_coeff/4.0 - h_b
    vertex_weights[4, bond] = 0.5 * xy_coeff
    vertex_weights[5, bond] = vertex_weights[4, bond]

    return 0
