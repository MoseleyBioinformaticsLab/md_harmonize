cimport cython
import numpy as np
cimport numpy as np

def compare(element1, element2):
    if element1[0] == "!":
        return not compare(element1[1:], element2)
    elif element2[0] == "!":
        return not compare(element1, element2[1:])
    elif element1 == "*" and element2 != "H":
        return True
    elif element2 == "*" and element1 != "H":
        return True
    return element1 == element2

@cython.boundscheck(False)
def find_mappings(np.ndarray[np.uint8_t, ndim=2] structure_matrix_1, np.ndarray[np.uint16_t, ndim=2] distance_matrix_1,
np.ndarray[np.uint8_t, ndim=2] structure_matrix_2, np.ndarray[np.uint16_t, ndim=2] distance_matrix_2, np.ndarray[np.uint8_t, ndim=2] mapping_matrix):

    cdef unsigned int num_rows = mapping_matrix.shape[0] - 1
    cdef unsigned int num_cols = mapping_matrix.shape[1] - 1
    cdef int row = 0
    cdef bint equal
    cdef np.ndarray[np.uint8_t, ndim=1] used = np.zeros(num_cols+1, dtype=np.uint8)
    cdef np.ndarray[np.long_t, ndim=1] mapping = np.array([-1 for _ in range(num_rows + 1)])
    cdef list isomorphs = []

    if not num_rows > num_cols:
        while row != -1:
            mapping[row] += 1
            while mapping[row] <= num_cols and (mapping_matrix[row, mapping[row]] == 0 or used[mapping[row]] == 1):
                mapping[row] += 1
            if not mapping[row] > num_cols:
                equal = True
                for j in range(row+1):
                    # the two atoms are connected and connected the same way in two substructures.
                    # Or the check the shortest distance between the two atoms in the more complex structure is shorter or equal than the simpler structure
                    if (structure_matrix_1[row, j] != 0 and structure_matrix_1[row, j] != structure_matrix_2[mapping[row], mapping[j]]) or (distance_matrix_1[row, j] < distance_matrix_2[mapping[row], mapping[j]]) :
                        equal = False
                        break
                if row == num_rows and equal:
                    isomorphs.append([x for x in mapping])
                elif row == 0 or equal:
                    used[mapping[row]] = 1
                    row += 1
            while row >= 0 and mapping[row] >= num_cols:
                mapping[row] = -1
                row -= 1
                used[mapping[row]] = 0
    return isomorphs

@cython.boundscheck(False)
def make_mapping_matrix(compound_1, compound_2, bint cyclic=False,  bint color_tuple=False, bint r_distance=False):

    cdef list atoms1
    cdef list atoms2
    cdef bint at_least_one
    cdef bint mappable
    cdef np.ndarray[np.uint8_t, ndim=2] mapping_matrix

    if compound_1.heavy_atoms and compound_2.heavy_atoms:
        atoms1 = compound_1.heavy_atoms
        atoms2 = compound_2.heavy_atoms
        mapping_matrix = np.zeros((len(atoms1), len(atoms2)), dtype=np.uint8)
        for i, atom_1 in enumerate(atoms1):
            at_least_one = False
            for j, atom_2 in enumerate(atoms2):
                mappable = compare(atom_1.default_symbol, atom_2.default_symbol)
                if cyclic and atom_1.in_cycle:
                    mappable = mappable and (atom_1.in_cycle == atom_2.in_cycle)
                if color_tuple:
                    mappable = mappable and all(x[0] <= x[1] for x in zip(atom_1.color_tuple, atom_2.color_tuple))
                if r_distance:
                    if atom_1.distance_to_r - 1 in atom_1.color_layers:
                        if atom_1.distance_to_r - 1 not in atom_2.color_layers or atom_1.color_layers[atom_1.distance_to_r - 1] != atom_2.color_layers[atom_1.distance_to_r - 1]:
                            mappable = False
                if mappable:
                    at_least_one = True
                mapping_matrix[i, j] = mappable
            if at_least_one is False:
                return None
    return mapping_matrix




