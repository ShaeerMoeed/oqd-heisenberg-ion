import numpy as np

def geometry_pbc(N, num_grid_pts):

    b = 0
    num_bonds = N * (num_grid_pts ** 2) # num_grid_pts^2 number of bonds between a single pair of rotors and N pairs of rotors
    geometry_table = np.zeros((num_bonds, 8)) # columns: b, i, j, i_x, i_y, j_x, j_y, r
    for i_x in range(N):
        for i_y in range(num_grid_pts):
            i = i_x*num_grid_pts + i_y

            j_x = i_x
            for j_y in range(i_y+1, num_grid_pts):
                j = j_x*num_grid_pts + j_y

                geometry_table[b, 0] = b
                geometry_table[b, 1] = i
                geometry_table[b, 2] = j
                geometry_table[b, 3] = i_x
                geometry_table[b, 4] = i_y
                geometry_table[b, 5] = j_x
                geometry_table[b, 6] = j_y
                geometry_table[b, 7] = np.sqrt((j_x - i_x)**2 + (j_y - i_y)**2)
                  
                b += 1

            j_x = (i_x + 1) % N
            for j_y in range(i_y, num_grid_pts):
                j = j_x*num_grid_pts + j_y

                geometry_table[b, 0] = b
                geometry_table[b, 1] = i
                geometry_table[b, 2] = j
                geometry_table[b, 3] = i_x
                geometry_table[b, 4] = i_y
                geometry_table[b, 5] = j_x
                geometry_table[b, 6] = j_y
                geometry_table[b, 7] = np.sqrt((j_x - i_x)**2 + (j_y - i_y)**2)

                b += 1

    return geometry_table