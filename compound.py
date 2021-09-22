#!/usr/bin/python3


"""
MDH.compound
~~~~~~~~~~~~~~~~~~~~~~~~~

This module provides the :class:`~MDH.compound.Atom` class, the :class:`MDH.compound.Bond` class,
and the :class:`~MDH.compound.Compound` class to construct compound entity.

"""

import collections
from math import atan, sin, cos
import itertools
import numpy
import BASS_Cc_3 as BASS_1
import heapq


class Atom:

	""" Atom class describes the :class:`~MDH.compound.Atom` entity in the compound. """

    def __init__(self, x, y, z, atom_symbol, mass_difference, charge, atom_stereo_parity, hydrogen_count, stereo_care_box, 
    	valence, h0designator, atom_atom_mapping_number, inversion_retention_flag, atom_number, exact_change_flag, kat=""):
    	"""Atom initializer.

        :param str x: the atom x coordinate.
        :param str y: the atom y coordinate.
        :param str z: the atom z coordinate.
        :param str atom_symbol: atom_symbol.
        :param str mass_difference: difference from mass in periodic table.
        :param str charge: charge.
        :param str atom_stereo_parity: atom stereo parity.
        :param str hydrogen_count: hydrogen_count.
        :param str stereo_care_box: stereo_care_box.
        :param str valence: valence.
        :param str h0designator: h0designator.
        :param str atom_atom_mapping_number: atom_atom_mapping_number.
        :param str inversion_retention_flag: inversion_retention_flag.
        :param str atom_number: atom_number.
        :param str exact_change_flag: exact_change_flag.
        :param str default_symbol: default_symbol.
        :param str kat: KEGG atom type.
    	"""
        self.x = float(x.strip())
        self.y = float(y.strip())
        self.z = float(z.strip())
        self.atom_symbol = atom_symbol.strip()
        self.mass_difference = mass_difference.strip()
        self.charge = int(charge.strip())
        self.atom_stereo_parity = atom_stereo_parity.strip()
        self.hydrogen_count = hydrogen_count.strip()
        self.stereo_care_box = stereo_care_box.strip()
        self.valence = valence.strip()
        self.h0designator = h0designator.strip()
        self.atom_atom_mapping_number = atom_atom_mapping_number.strip()
        self.inversion_retention_flag = inversion_retention_flag.strip()
        self.atom_number = atom_number
        self.exact_change_flag = exact_change_flag.strip()
        self.neighbors = []
        self.color_0 = ""
        self.color = ""
        self.color_tuple = tuple()
        self.color_layers = {}
        self.in_cycle = False
        self.bond_counts = 0
        self.group_id = 0
        self.double_bond_counts = 0
        self.distance_to_R = 0
        self.kat = kat
        self.atom_mappings = []

    @property
    def default_symbol(self):

        if self.is_R():
            return "R"
        return self.atom_symbol

    def update_symbol(self, symbol):
    	"""
        To update the symbol of the atom.
        :param str symbol: the updated atom symbol.

        :return: the updated atom_symbol.
    	"""
        self.atom_symbol = symbol
        return self.atom_symbol

    def remove_neighbors(self, neighbors):
    	"""
        To add neighbors to the atom.
        :param list neighbors: the index of neighbors that will be removed from this atom.

        :return: the updated list of neighbors of the atom.
        """
        for atom_index in neighbors:
            if atom_index in self.neighbors:
                self.neighbors.remove(atom_index)
        return self.neighbors

    def add_neighbors(self, neighbors):
        """
        To add neighbors to the atom.
        :param list neighbors: the index of neighbors that will be added to this atom.

        :return: the updated list of neighbors of the atom.
        """
        for atom_index in neighbors:
            if atom_index not in self.neighbors:
                self.neighbors.append(atom_index)
        return self.neighbors

    def update_stereochemistry(self, stereo):
    	"""
        To update the stereochemistry of the atom.
        :param str stereo: the updated atom stereochemistry.

        :return: the updated stereochemistry.
    	"""
        self.atom_stereo_parity = stereo
        return self.atom_stereo_parity
    
    def color_atom(self, isotope_resolved=False, charge=False, atom_stereo=False):
    	"""
        To generate the zero layer atom color.
        :param isotope_resolved: If true, add isotope information when constructing colors.
        :type isotope_resolved: :py:obj:`True` or :py:obj:`False`.
        :param charge: If true, add charge information when constructing colors.
        :type charge: :py:obj:`True` or :py:obj:`False`.
        :param atom_stereo: If true, add atom stereochemistry information when constructing colors.
        :type atom_stereo: :py:obj:`True` or :py:obj:`False`.

	    :return: the zero layer atom color.
    	"""
    	self.color_0 = ""
        self.color_0 += self.default_symbol
        if atom_stereo:
            self.color_0 += self.atom_stereo_parity
        if charge:
            self.color_0 += str(self.charge)
        if isotope_resolved:
            self.color_0 += self.mass_difference
        self.color = self.color_0
    	return self.color_0

    def reset_color(self):
    	"""
        Set the color of the atom to empty.

        :return:
    	"""
    	self.color_0 = ""
    	self.color = ""
    	self.color_layers = {}
    
    def is_R(self):
    	"""
        To determine if the atom is an R group.

        :return: boolean.
    	"""
        if "A" in self.default_symbol or "R" in self.default_symbol or "*" in self.default_symbol:
            if self.default_symbol not in not_R_groups:
                return True
        return False


class Bond:

	""" Bond class describes the :class:`~MDH.compound.Bond` entity in the compound. """

    def __init__(self, first_atom_number, second_atom_number, bond_type, bond_stereo, bond_topology, reacting_center_status):
    	"""Bond initializer.

        :param str first_atom_number: the index of the first atom forming this bond.
        :param str second_atom_number: the index of the second atom forming this bond.
        :param str bond_type: the bond type. (1 = Single, 2 = Double, 3 = Triple, 4 = Aromatic, 5 = Single or Double,
        6 = Single or Aromatic, 7 = double or Aromatic 8 = Any)
        :param str bond_stereo: the bond stereo. (Single bonds: 0 = not stereo, 1 = Up, 4 = Either, 6 = Down; Double bonds: determined by x, y, z coordinates)
        :param str bond_topology: bond topology. (O = Either, 1 = Ring, 2 = Chain)
        :param str reacting_center_status: reacting center status.
    	"""
        self.first_atom_number = int(first_atom_number) - 1
        self.second_atom_number = int(second_atom_number) - 1
        self.bond_type = bond_type.strip()
        self.bond_stereo = bond_stereo.strip()
        self.bond_topology = bond_topology.strip()
        self.reacting_center_status = reacting_center_status.strip()
 
    def update_bond_type(self, bond_type):
    	"""
        To update the bond type of the atom.
        :param str bond_type: the updated bond type.

        :return: the updated bond type.
    	"""
        self.bond_type = str(bond_type)
        return self.bond_type

    def update_stereochemistry(self, stereo):
    	"""
        To update the stereochemistry of the bond.
        :param str stereo: the updated atom stereochemistry.

        :return: the updated stereochemistry.
    	"""
        self.bond_stereo = str(stereo)
        return self.bond_stereo


class Compound:

	""" Compound class describes the :class:`~MDH.compound.Compound` entity. """

    def __init__(self, compund_name, atoms, bonds):
    	"""Compound initializer.
			
        :param str compund_name: the compund name.
        :param list atoms: a list of :class:`~MDH.compound.Atom` enties in the compound.
        :param list bonds: a list of :class:`~MDH.compound.Bond` enties in the compound.
    	"""    
    	self.compund_name = compund_name
        self.atoms = atoms
        self.bonds = bonds
        self.bond_lookup = {}
        self.has_cycle = False
        
        for bond in self.bonds:
            if bond.bond_type != "8":
                first_atom, second_atom = self.atoms[bond.first_atom_number], self.atoms[bond.second_atom_number]
                self.bond_lookup[(bond.first_atom_number, bond.second_atom_number)] = bond
                self.bond_lookup[(bond.second_atom_number, bond.first_atom_number)] = bond
                first_atom.neighbors.append(bond.second_atom_number)
                first_atom.bond_counts += int(bond.bond_type)
                second_atom.neighbors.append(bond.first_atom_number)
                second_atom.bond_counts += int(bond.bond_type)     
        
        self.find_cycles()
        self.calculate_distance_to_R_groups()

    @property
    def formula(self):
    	"""
        To construct the formula of this compound (only consider heavy atoms).

        :return: string.
    	"""
        counter = collections.Counter()
        for atom in self.atoms:
            if atom.default_symbol == "H":
                continue
            elif atom.is_R():
                counter["R"] += 1
            else:
                counter[atom.default_symbol] += 1
        return "".join([ char + str(counter[char]) for char in sorted(counter.keys()) ])
    
    @property
    def R_groups(self):
    	"""
        To get all the R groups in the compound.

        :return: the list of index of all the R groups.
    	"""
        return [ index for index, atom in enumerate(self.atoms) if atom.is_R() ]

    def contains_R_groups(self):
    	"""
        To check if the compound contains R group(s).

        :return: boolean.
    	"""
        return self.R_groups != []

    def has_isolated_atoms(self):
    	"""
        To check if the compound has atoms that have no connections to other atoms.

        :return: boolean
    	"""
        for atom in self.atoms:
            if atom.bond_counts == 0:
                return True
        return False

    @property
    def metal_index(self):
        """
        To get the indexes of metals without connections to other atoms in the compound.

        :return: the list of index of isolated metals.
        """
        return [index for index, atom in enumerate(self.atoms) if atom.default_symbol in metal_symbols]
    
    @property
    def H_index(self):
    	"""
        To get the index of H(s) in the compound.

        :return: the list of index of Hs.
    	"""
        return [index for index, atom in enumerate(self.atoms) if atom.default_symbol == "H"]

    @property
    def heavy_atoms(self):
        """
        To get all the heavy atoms in the compound.
        :param self:
        :return: the list of heavy atoms in the compound.
        """

        return [atom for atom in self.atoms if atom.default_symbol != "H"]
    
    @property
    def index_of_heavy_atoms(self):
    	"""
        To map the atom index to heavy atoms.

        :return: the dictionay of atom index to atom entity.
    	"""
        return { atom.atom_number: i for i, atom in enumerate(self.heavy_atoms) }

       
    def color_groups(self, excluded=None):
    	"""
        To update the compound color groups after coloring.

        :return: the dictionary of atom color name with the corresponding atom index list.
    	"""
    	color_groups = collections.defaultdict(list)
        for atom in self.atoms:
            if atom.atom_number not in excluded:
                color_groups[atom.color].append(atom.atom_number)
        return color_groups

    def detect_abnormal_atom(self):
    	"""
        Find the atoms with invalid bond counts.

        :return: the list of abnomal atom indexes.
    	"""
        abnormal_atoms = collections.defaultdict(list)
        for index, atom in enumerate(self.atoms):
            if atom.default_symbol in standard_bond_counts:
                bond_counts = atom.bond_counts
                bond_counts -= atom.charge
                if type(standard_bond_counts[atom.default_symbol]) is list:
                    if bond_counts not in standard_bond_counts[atom.default_symbol]:
                        abnormal_atoms[atom.default_symbol].append(index)
                else:
                    if bond_counts > standard_bond_counts[atom.default_symbol]:
                        abnormal_atoms[atom.default_symbol].append(index)
        return abnormal_atoms

    def curate_invalid_N(self):
    	"""
        Curate the invalid N atoms.

        :return:
    	"""
        abnormal_atoms = self.detect_abnormal_atom()
        for atom_index in abnormal_atoms["N"]:
            atom = self.atoms[atom_index]
            atom.charge += 1
            for neighbor in atom.neighbors:
                if self.atoms[neighbor].default_symbol == "O" and self.bond_lookup[(atom.atom_number, neighbor)].bond_type == "2":
                    self.atoms[neighbor].charge -= 1
                    self.bond_lookup[(atom.atom_number, neighbor)].update_stereochemistry("1")

    def metals(self):
    	"""
        To get all the metals in the compound.

        :return: the dictionary of metal and the corresponding list of index.
    	"""
        metals = collections.defaultdict(list)
        for index, atom in enumerate(self.atoms):
            if atom.default_symbol in metal_symbols:
                metals[atom.default_symbol].append(index)
        return metals

 #    def update_double_bond_stereo(self, bonds):
	# """
	# To update the stereochemistry of the double bonds.
	# :param list bonds: the list of bonds whose bond stereo needs update. 

	# :return: 
	# """
	# for bond in bonds:
	#     first_atom = bond.first_atom_number
	#     second_atom = bond.second_atom_number
	#     if self.bond_lookup[(first_atom, second_atom)].bond_type == "2":
	#         self.bond_lookup[(first_atom, second_atom)].update_stereochemistry(bond.bond_stereo)

    def update_aromatic_bond_type(self, cycles):
        """
        Update the aromatic bond types.
        Two cases: 1) change the bond in the aromatic ring to aromatic bond (bond type = 4)
               2) change the double bond connecting to the aromatic ring to single bond.
        :param list cycles: the list of cycls of aromatic atom indexes.

        :return:
        """
        atom_in_cycle = [atom for cycle in cycles for atom in cycle]
        for cycle in cycles:
            aromatic_bonds = self.extract_aromatic_bonds(cycle)
            for bond in aromatic_bonds:
                bond.update_bond_type("4")
        bond_out_of_cycle = self.extract_double_bond_connecting_cycle(atom_in_cycle)
        for bond in bond_out_of_cycle:
            bond.update_bond_type("1")

    def extract_double_bond_connecting_cycle(self, atom_in_cycle):
        """
        Extract the double bonds connecting to the aromatic cycles.
        :param list atom_in_cycle: the list of indexes of atoms in the aromatic cycles.

        :return: the list of double bonds connecting to the aromatic cycles.
        """
        double_bond_connecting_cycle = []
        for atom_index in atom_in_cycle:
            if self.atoms[atom_index].default_symbol == "C":
                for neighbor_index in self.atoms[atom_index].neighbors:
                    if neighbor_index not in atom_in_cycle:
                        if self.bond_lookup[(atom_index, neighbor_index)].bond_type == "2":
                            double_bond_connecting_cycle.append(self.bond_lookup[(atom_index, neighbor_index)])
        return double_bond_connecting_cycle

    def extract_aromatic_bonds(self, cycle):
    	"""
        Extract the aromatic bonds based on the atom indexes in the cycle.
        :param list cycle: the

        :return: the list of aromatic bonds.
    	"""
        aromatic_bonds = []
        all_pairs = list(itertools.combinations(cycle, 2))
        visited = set()
        for pair in all_pairs:
            if (pair[0], pair[1]) not in visited and (pair[0], pair[1]) in self.bond_lookup:
                aromatic_bonds.append(self.bond_lookup[(pair[0], pair[1])])
                visited.add((pair[0], pair[1]))
                visited.add((pair[1], pair[1]))
        return aromatic_bonds

    def connected_components(self):
    	"""
        Detect the connected components in the compound structure. (using the breadth first search)

        :return: the dictionary of the connected components.
    	"""
        visited = set()
        group_id = 0
        components = collections.defaultdict(list)
        for index, atom in enumerate(self.atoms):
            if index not in visited:
            	ques = collections.deque(index)
            	while ques:
            		cur = ques.popleft()
            		self.atoms[cur].group_id = group_id
            		visited.add(cur)
            		components[group_id].append(cur)
            		for neighbor in self.atoms[cur].neighbors:
            			if neighbor not in visited:
            				ques.append(neighbor)
                group_id += 1
        return components

    def calculate_distance_to_R_groups(self):
        """
        To caluclate the distance of each atom to its nearest R group (using the dijkstra's algorithm).

        :return:
        """
        distance_matrix = [len(self.heavy_atoms)] * len(self.heavy_atoms) 

        if self.R_groups:
            for r_index in self.R_groups:
                ques = [(0, r_index)]
                seen = {}
                while ques:
                    dist, atom_index = heapq.heappop(ques)
                    if atom_index in seen and seen[atom_index] <= dist:
                        continue
                    seen[atom_index] = dist
                    distance_matrix[self.index_of_heavy_atoms[atom_index]] = min(distance_matrix[self.index_of_heavy_atoms[atom_index]], dist)
                    for neighbor in self.atoms[atom_index].neighbors:
                        if self.atoms[neighbor].default_symbol != "H":
                            if neighbor not in seen or seen[neighbor] > dist+ 1:
                                if distance_matrix[self.index_of_heavy_atoms[neighbor]] > dist+1:
                                    heapq.heappush(ques, (dist+1, neighbor))
        
            for i, dist in enumerate(distance_matrix):
                self.heavy_atoms[i].distance_to_R = dist

    def find_cycles(self, short_circuit=False, cutoff=40):
    	"""
        Find the cycles in the compound.

        :return:
    	"""
        atoms, not_cyclic, cyclic, all_cycles = self.atoms, [], [], []

        def prune():
            while 1:
                terminate = True
                for prune_atom in atoms:
                    if prune_atom.atom_number not in not_cyclic:
                        # for atom with cycle, it at least has 2 neighbors to form a cycle.
                        if len([neighbor for neighbor in prune_atom.neighbors if neighbor not in not_cyclic]) < 2:
                            not_cyclic.append(prune_atom.atom_number)
                            terminate = False
                if terminate:
                    break

        def search(start_index):
            # the paths is a list of list, which means each path is a list.
            paths, cycles = collections.deque([[start_index]]), []
            while paths:
                indices = paths.popleft()
                for neighbor in atoms[indices[-1]].neighbors:
                    if neighbor not in not_cyclic:
                        if neighbor not in indices and len(indices) < cutoff:
                            paths.append(list(indices) + [atoms[neighbor].atom_number])
                        elif neighbor == start_index and len(indices) != 2:
                            cycles.append(list(indices) + [atoms[neighbor].atom_number])
                            if short_circuit:
                                return cycles
            if cycles:
                return cycles
        prune()
        for atom in atoms:
            start_index = atom.atom_number
            if start_index not in not_cyclic:
                if not short_circuit or start_index not in cyclic:
                    search_result = search(start_index)
                    if search_result is None:
                        not_cyclic.append(start_index)
                        prune()
                    else:
                        all_cycles.append(search_result)
                        if short_circuit:
                            cyclic += [item for sublist in search_result for item in sublist]

        for cycle_list in all_cycles:
            for cycle in cycle_list:
                for index in cycle:
                    atoms[index].in_cycle = True
        
        if all_cycles:
            self.has_cycle = True

    @property
    def resonance_targeted_structure_matrix(self):
    	"""
        To construct graph structural matrix of this compound without distinguishing single and double bonds.
        This is used to find compounds with different resonance structures.

        :return: the matrix of this compound.
    	"""
        resonance_targeted_structure_matrix = numpy.zeros((len(self.heavy_atoms), len(self.heavy_atoms)), dtype=numpy.uint8)
        for bond in self.bonds:
            atom_1, atom_2 = bond.first_atom_number, bond.second_atom_number
            if atom_1 in self.index_of_heavy_atoms and atom_2 in self.index_of_heavy_atoms:
                resonance_targeted_structure_matrix[self.index_of_heavy_atoms[atom_1]][self.index_of_heavy_atoms[atom_2]] = 
                resonance_targeted_structure_matrix[self.index_of_heavy_atoms[atom_2]][self.index_of_heavy_atoms[atom_1]] = int(bond.bond_type) if int(bond.bond_type) != 2 else 1
        return resonance_targeted_structure_matrix

    @property
    def structure_matrix(self):
    	"""
        To construct graph structural matrix of this compound.
        matrix[i][j] = 0 suggests the two atoms are not connected directly.

        :return: the structure matrix of this compound.
    	"""
        matrix = numpy.zeros((len(self.heavy_atoms), len(self.heavy_atoms)), dtype=numpy.uint8)
        for bond in self.bonds:
            atom_1, atom_2 = bond.first_atom_number, bond.second_atom_number
            if atom_1 in self.index_of_heavy_atoms and atom_2 in self.index_of_heavy_atoms:
                matrix[self.index_of_heavy_atoms[atom_1]][self.index_of_heavy_atoms[atom_2]] = 
                matrix[self.index_of_heavy_atoms[atom_2]][self.index_of_heavy_atoms[atom_1]] = int(bond.bond_type)
        return matrix

    @property
    def distance_matrix(self):
    	"""
        To construct the distance matrix of the compound. (using the Floyd Warshall Algorithm)
        distance[i][j] suggests the distance between atom i and j.

        :return: the distance matrix of the compound.
    	"""
        if self.heavy_atoms:
            distance_matrix = numpy.ones((len(self.heavy_atoms), len(self.heavy_atoms)),dtype=numpy.uint16 ) * len(self.heavy_atoms)
            for bond in self.bonds:
                atom_1, atom_2 = bond.first_atom_number, bond.second_atom_number
                if atom_1 in self.index_of_heavy_atoms and atom_2 in self.index_of_heavy_atoms:
                    distance_matrix[self.index_of_heavy_atoms[atom_1]][self.index_of_heavy_atoms[atom_2]] = 
                    distance_matrix[self.index_of_heavy_atoms[atom_2]][self.index_of_heavy_atoms[atom_1]] = 1

            n = len(self.heavy_atoms)
            for k in range(n):
                for i in range(n):
                    for j in range(n):
                        if distance_matrix[i][j] > distance_matrix[i][k] + distance_matrix[k][j]:
                            distance_matrix[i][j] = distance_matrix[i][k] + distance_matrix[k][j]
            return distance_matrix
        return None
    
    def update_color_tuple(self, resonance=False):
    	"""
        To update the color tuple of the atoms in the compound. This color tuple includes information of its neighboring atoms and bonds.
        :param resonance: to find the resonant structure. This will ignore the difference between single and double bonds.
        :type resonance: :py:obj:`True` or :py:obj:`False`.

        :return:
    	"""
        for atom in self.heavy_atoms:
            elements = collections.Counter()
            bond_types = collections.Counter()
            for neighbor_index in atom.neighbors:
                neighbor = self.atoms[neighbor_index]
                if neighbor.default_symbol != "H":
                    elements[neighbor.default_symbol] += 1
                    bond = self.bond_lookup[(atom.atom_number, neighbor_index)]
                    bond_types[bond.bond_type] += 1
            if resonance:
            	atom.color_tuple = (elements["C"], elements["N"], elements["S"], elements["O"], 
            		bond_types["1"] + bond_types["2"], bond_types["4"])
            else:
            	atom.color_tuple = (elements["C"], elements["N"], elements["S"], elements["O"], 
            		bond_types["1"], bond_types["2"], bond_types["4"])
    
 #    def update_mappings_atoms(self, the_other, mappings):
 #    	"""
	# Update the mappings atoms between two compounds.
	# :param the_other: the mappings compound entity.
	# :type the_other: :class:`~MDH.compound.compound`
	# :param list mappings: the mappings atom indexes between the two compounds.

	# :return:
 #    	"""
 #        for atom in self.atoms:
 #            atom.atom_mappings = []
 #        for atom in the_other.atoms:
 #            atom.atom_mappings = []

 #        if mappings:
 #            for atom_1 in mappings:
 #                for atom_2 in mappings[atom_1]:
 #                    self.atoms[atom_1].atom_mappings.append(atom_2)
 #                    the_other.atoms[atom_2].atom_mappings.append(atom_1)

   

    def find_mappings(self, the_other, resonance=False, R=False, default_mappings=None):
    	"""
        Find the atom mappings between two compounds.
        :param the_other: the mappings compound entity.
        :type the_other: :class:`~MDH.compound.Compound`
        :param resonance: to find the resonant structure. This will ignore the difference between single and double bonds.
        :type resonance: :py:obj:`True` or :py:obj:`False`.
        :param R: to take account of positon of R groups.
        :type R: :py:obj:`True` or :py:obj:`False`.
        :param list default_mappings: the default mappings between the two compounds.

        :return: the list of atom mappings in the heavy atom order.
    	"""
        mappping_matrix = BASS_1.make_mapping_matrix_R if R else BASS_1.make_mapping_matrix 
        find_mappings = BASS_1.find_mappings
        
        self.update_color_tuple(resonance=resonance)
        the_other.update_color_tuple(resonance=resonance)
        mappings = []
        # self.update_mappings_atoms(the_other, default_mappings)

        mmat = mappping_matrix(the_other, self, True, True)
        if mmat is not None:
            mappings = find_mappings(the_other.structure_matrix, the_other.distance_matrix, 
            	self.structure_matrix, self.distance_matrix, mmat)
        return mappings

    def generate_one_to_one_mappings(self, the_other, mappings):
    	"""
        Generate the one to one atom mappings between two compounds.
        :param the_other: the mappings compound entity.
        :type the_other: :class:`~MDH.compound.Compound`
        :param list mappings: the mappings of heavy atom indexes.

        :return: the list of dictionary of one to one atom mappings in the original order.
    	"""
    	one_to_one_mappings = []
        for sub in mappings:
            cur_mappings = {}
            for from_index, to_index in enumerate(sub):
                cur_mappings[self.heavy_atoms[to_index].atom_number] = the_other.heavy_atoms[from_index].atom_number
            one_to_one_mappings.append(cur_mappings)
        return one_to_one_mappings

    def validate_resonance_mappings(self, the_other, R=False):
    	"""
        Check if the resonant mappings are valid between the two compound structures. If the mapped atoms don't share the same
        local coloring identifier, we check if the difference is caused by the position of double bonds.
        Find the three atoms involved in the resonant structure and check if one of the atom is not C.
                N (a)            N (a)
		        / \\            // \
            (b) C   N (c)    (b) C   N (c)

        :param the_other: the mappings compound entity.
        :type the_other: :class:`~MDH.compound.Compound`
        :param R: to take account of positon of R groups.
        :type R: :py:obj:`True` or :py:obj:`False`.

	    :return: the list of valid atom mappings between the two compound structures.
    	"""
    	mappings = self.find_mappings(self, the_other, resonance=True, R=False)
        one_to_one_mappings = self.generate_one_to_one_mappings(the_other, mappings)

        valid_mappings = []
        for cur_mappings in one_to_one_mappings:
            flag = True
            reversed_mappings = { cur_mappings[key] : key for key in cur_mappings }
            for from_index in sorted(cur_mappings.keys()):
                to_index = cur_mappings[from_index]

                "Detect if the two atoms have the same local binding environment by comparing the first layer color identifier."
                if 1 in self.atoms[from_index].color_layers and 1 in the_other.atoms[to_index].color_layers and 
                self.atoms[from_index].color_layers[1] == the_other.atoms[to_index].color_layers[1]:
                    continue
                
                three = {from_index}
                from_dbond_atom_index = self.find_double_bond_linked_atom(from_index)
                to_dbond_atom_index = the_other.find_double_bond_linked_atom(to_index)
                
                "cannot find directly linked double bond."
                if from_dbond_atom_index == -1 and to_dbond_atom_index == -1:
                    flag = False
                    break

                "Need to conduct second search"
                "atom b in the picture"
                elif from_dbond_atom_index == -1:
                    reversed_from_dbond_atom_index = reversed_mappings[to_dbond_atom_index] # here we can find the atom a
                    three.add(reversed_from_dbond_atom_index)
                    from_next_dbond_atom_index = self.find_double_bond_linked_atom(reversed_from_dbond_atom_index) 
                    # based on atom a we can find atom c
                    if from_next_dbond_atom_index != -1:
                        three.add(from_next_dbond_atom_index)

                "atom c in the picture"
                elif to_dbond_atom_index == -1 and from_dbond_atom_index in cur_mappings:
                    three.add(from_dbond_atom_index) # here we find atom a
                    reversed_to_dbond_atom_index = cur_mappings[from_dbond_atom_index]
                    # based on atom a find atom b in the other compound.
                    to_next_dbond_atom_index = the_other.find_double_bond_linked_atom(reversed_to_dbond_atom_index)
                    if to_next_dbond_atom_index != -1:
                        three.add(reversed_mappings[to_next_dbond_atom_index])

                "atom a in the picture"       	
                "find the other two atoms directly."
                elif from_dbond_atom_index != -1 and to_dbond_atom_index != -1:
                    three.add(from_dbond_atom_index)
                    three.add(reversed_mappings[to_dbond_atom_index])
                
                if len(three) < 3:
                    flag = False
                    break
                
                k = [self.atoms[i].default_symbol for i in three].count("C")
                if k == 3:
                    flag = False
                    break
            if flag:
                valid_mappings.append(cur_mappings)
        return valid_mappings
    
    def find_double_bond_linked_atom(self, index):
    	"""
        Find the index of the atom that is doubly linked to the target atom[i].
        param int index: the index of the target atom.

        :return: the index of the doubly linked atom.
    	"""
        for neighbor_index in self.atoms[index].neighbors:
            if neighbor_index in self.index_of_heavy_atoms and self.bond_lookup[(index, neighbor_index)].bond_type == "2":
                return neighbor_index
        return -1

    def detect_linear_circular_interchange(self, the_other):
    	"""
        Detect if the two compounds have linear and circular interchangeable representations. This mainly targets at sugar like structure.
        We check how many atoms don't share the same local (first layer) identifiers and then detect if the difference is casued by the
        ring structure.
        :param the_other: the mappings compound entity.
        :type the_other: :class:`~MDH.compound.Compound`

        :return: boolean.
    	"""
    	# 1 in color_layers indicate this atom is not H or metals.
        one_atom_colors = [atom.color_layers[1] if 1 in atom.color_layers else atom.default_symbol for atom in self.heavy_atoms]
        the_other_atom_colors = [atom.color_layers[1] if 1 in atom.color_layers else atom.default_symbol for atom in self.heavy_atoms]
        atom_color_intersection = list((collections.Counter(one_atom_colors) & collections.Counter(the_other_atom_colors)).elements())
        if 2 <= len(one_atom_colors) - len(atom_color_intersection) <= 3:
            one_num_of_atoms_in_cycle = [ atom.in_cycle for atom in self.heavy_atoms].count(True)
            the_other_num_of_atoms_in_cycle = [ atom.in_cycle for atom in the_other.heavy_atoms].count(True)
            if one_num_of_atoms_in_cycle != the_other_num_of_atoms_in_cycle:
                return True
        return False

    def define_bond_stereochemistry(self):
    	"""
        Define the stereochemistry of double bonds in the compounds.

        :return:
    	"""
        for bond in self.bonds:
            if bond.bond_type == "2":
            	bond.update_stereochemistry(self.calculate_bond_stereochemistry(bond))
    
    def calculate_bond_stereochemistry(self, bond):
    	"""
        Calculate the stereochemisty of double bond based on the geometric properties. The line of double bond divivdes the plane into two parts. For the atom forming the
        doble bond, it normally has two branches. If the two branches are not the same, we can them heavy side and light side (heavy side containing atoms with heavier atomic weights).
        We determine the bond stereochemisty by checking if the two heavy sides lie on the same part of the divided plane.

           H     L     H     H
            \___/	    \___/
             ___         ___
            /	\       /   \
           L 	 H     L     L
            trans        cis
        :param bond: the bond entity.
        :type bond: :class:`~MDH.compound.Bond`

        :return: the calculated bond stereochemistry.
    	"""
        vertical = False

        first_atom = self.atoms[bond.first_atom_number]
        second_atom = self.atoms[bond.second_atom_number]

        if first_atom.x != second_atom.x:
            slope = (first_atom.y-second_atom.y)/(first_atom.x-second_atom.x)
            b = first_atom.y - slope * first_atom.x 
        else:
            vertical = True

        first_atom_neighor_index = [neighbor for neighbor in first_atom.neighbors if neighbor != second_atom.atom_number]
        second_atom_neighor_index = [neighbor for neighbor in second_atom.neighbors if neighbor != first_atom.atom_number]
        
        if not first_atom_neighor_index or not second_atom_neighor_index:
            return 0
        
        first_heavy_side, first_light_side = self.rank_sides(first_atom_neighor_index, first_atom)
        second_heavy_side, second_light_side = self.rank_sides(second_atom_neighor_index, second_atom)

        if (first_light_side and first_light_side == first_heavy_side) or 
        (second_light_side and second_light_side == second_heavy_side):
        	return 0
       
        if vertical:
            if (self.atoms[first_heavy_side].x - first_atom.x) * (self.atoms[second_heavy_side].x - second_atom.x) > 0:
                return 1
            else:
                return -1
        else:
            if (self.calculate_y_coordinate(slope, b, self.atoms[first_heavy_side]) - self.atoms[first_heavy_side].y) 
            * (self.calculate_y_coordinate(slope, b, self.atoms[second_heavy_side]) - self.atoms[second_heavy_side].y) > 0:
                return 1
            else:
                return -1

    def calculate_y_coordinate(self, slope, b, atom):
    	"""
        Calculate the y coordinate of the atom based on linear function. y = slope * x + b
        :param float slope: the slope of the targeted line.
        :param float b: the intercept of the targeted line.

        :return: the calculated y coordinate.
    	"""
        return atom.x * slope + b

    def collect_atomic_weights_of_neighbors(self, neighbors):
    	"""
        To collect the atomic weights of the next layer's neighbors.
        :param list neighbors: the list of this layer's atom with its parent atom. The parent is used to avoid recalculating the former layer.

        :return: the list of atomic weights of next layer's neighbors
    	"""
        neighbor_atomic_weights = [atomic_weights[self.atoms[index].default_symbol] for index in neighbors]
        neighbor_atomic_weights.sort(reverse=True)
        return neighbor_atomic_weights

    def rank_sides(self, neighbors, atom_forming_double_bond):
    	"""
        To determine the rank of the two branches connecting the atom forming the double bond. This is based on comparison of the atomic weights of the two branches.
        Breadth first algorithm.
        :param list neighbors: the list of atom indexed of the atoms connecting the atom forming the double bond.
        :param atom_forming_double_bond: the bond entity.
	    :type atom_forming_double_bond: :class:`~MDH.compound.Atom`
    	
        :return: the list ranked atoms. [heavy_side, light_side]
        """
    	if len(neighbors) < 2:
    		return self.atoms[neighbors[0]], None

    	one_atom, the_other_atom = neighbors[0], neighbors[1]

        if atomic_weights[self.atoms[one_atom].default_symbol] > atomic_weights[self.atoms[the_other_atom].default_symbol]:
            return [one_atom, the_other_atom]

        elif atomic_weights[self.atoms[one_atom].default_symbol] < atomic_weights[self.atoms[the_other_atom].default_symbol]:
            return [the_other_atom, one_atom]

        else:
            one_neighbors = [one_atom]
            one_visited = {atom_forming_double_bond.atom_number}
            the_other_neighbors = [the_other_atom]
            the_other_visited = {atom_forming_double_bond.atom_number}
            one_neighbor_atomic_weight_list = self.collect_atomic_weights_of_neighbors(one_neighbors)
            the_other_neighbor_atomic_weight_list = self.collect_atomic_weights_of_neighbors(the_other_neighbors)

            while one_neighbor_atomic_weight_list == the_other_neighbor_atomic_weight_list and one_neighbor_atomic_weight_list:
                one_neighbors = self.get_next_layer_neighbors(one_neighbors, one_visited)
                one_neighbor_atomic_weight_list = self.collect_atomic_weights_of_neighbors(one_neighbors)
                the_other_neighbors = self.get_next_layer_neighbors(the_other_neighbors, the_other_visited)
                the_other_neighbor_atomic_weight_list = self.collect_atomic_weights_of_neighbors(the_other_neighbors)

            if tuple(one_neighbor_atomic_weight_list) > tuple(the_other_neighbor_atomic_weight_list):
                return [one_atom, the_other_atom]
            elif tuple(one_neighbor_atomic_weight_list) < tuple(the_other_neighbor_atomic_weight_list):
                return [the_other_atom, one_atom]
            else:
                return [one_atom, one_atom]

    def get_next_layer_neighbors(self, cur_layer_neighbors, visited):
    	"""
        To get the next layer's neighbors.
        :param list cur_layer_neighbors: the list of atom index of neighoring atoms at the current layer.
        :param set visited: the set of atom indexes that have already been visited.

        :return: the next layer's neighbors.
    	"""
    	next_layer_neighbors = []
    	for index in cur_layer_neighbors:
    		for next_neighbor in self.atoms[index].neighbors:
    			if next_neighbor not in visited:
    				next_layer_neighbors.append(next_neighbor)
    				visited.add(next_neighbor)
    	return next_layer_neighbors

    def color_compound(self, R=True, bond_stereo=True, atom_stereo=True, resonance=True):
    	"""
        Color the compond.
        :param R:  If true, add R groups in the coloring.
        :type R: :py:obj:`True` or :py:obj:`False`.
        :param bond_stereo:  If true, add stereo information to bonds when constructing colors.
        :type bond_stereo: :py:obj:`True` or :py:obj:`False`.
        :param atom_stereo: If true, add atom stereo information when constructing colors.
        :type atom_stereo: :py:obj:`True` or :py:obj:`False`.
        :param resonance: If true, detect resonant compound pairs without disginguish between double bonds and single bonds.
        :type resonance: :py:obj:`True` or :py:obj:`False`.

        :return:
    	"""
        self.first_round_color(R=R, bond_stereo=bond_stereo, atom_stereo=atom_stereo, resonance=resonance)
        self.curate_invalid_symmetric_atoms(bond_stereo=bond_stereo, resonance=resonance)
        self.color_metals(bond_stereo=bond_stereo, resonance=resonance)
    
    def reset_color(self):
    	"""
        To set the color of atoms and compound to empty.

        :return:
    	"""
        for atom in self.atoms:
            atom.reset_color()

    def generate_atom_zero_layer_color(self, isotope_resolved=False, charge=False, atom_stereo=False):
    	"""
        To generate the zero layer color identifier for each atom. We don't consider H and metals here.
        :param isotope_resolved: If true, add isotope information when constructing colors.
        :type isotope_resolved: :py:obj:`True` or :py:obj:`False`.
        :param charge: If true, add charge information when constructing colors.
        :type charge: :py:obj:`True` or :py:obj:`False`.
        :param R: If true, add R group information when constructing colors.
        :type R: :py:obj:`True` or :py:obj:`False`.
        :param atom_stereo: If true, add atom stereochemistry information when constructing colors.
        :type atom_stereo: :py:obj:`True` or :py:obj:`False`.

        :return:
    	"""
    	for index, atom in enumerate(self.atoms):
            atom.color_atom(isotope_resolved=isotope_resolved, charge=charge, atom_stereo=atom_stereo)

    def generate_atom_color_with_neighbors(self, atom_index, excluded=None, zero_layer=True, resonance=False, bond_stereo=False):
    	"""
        To generate the atom color with its neighbors. We add this color name when we try to incorprate neighbors' information in naming.
        Here, we don't need to care about the atom stereo. It has been taken care of in generating color_0.
        Basic color formula: atom.color + [neighbor.color + bond.bond_type]
        :param zero_layer: If ture, we use the atom.color_0 else atom.color.
        :type zero_layer: :py:obj:`True` or :py:obj:`False`.
        :param resonance: If true, detect resonant compound pairs without disginguish between double bonds and single bonds.
        :type resonance: :py:obj:`True` or :py:obj:`False`.
        :param bond_stereo:  If true, add stereo information to bonds when constructing colors.
        :type bond_stereo: :py:obj:`True` or :py:obj:`False`.

        :return: the dictionary of atom index to its color name containing neighbors.
    	"""
    	atom_color_with_neighbors = collections.defaultdict(str)
    	for index in atom_index:
            atom_color = atom.color_0 if zero_layer else atom.color
            color_elements = collections.Counter()
           	for neighbor_index in atom.neighbors:
           		if neighobr_index not in excluded:
                    neighbor_color = self.atoms[neighbor_index].color_0 if zero_layer else self.atoms[neighbor_index].color
                    connecting_bond = self.bond_lookup[(atom.atom_number, neighbor_index)]
                    bond_type = connecting_bond.bond_type
                    if resonance and bond_type == "2":
                        bond_type = "1"
                    if bond_stereo:
                        component = "({0}.{1})".format(neighbor_color, bond_type+connecting_bond.bond_stereo)
                    else:
                        component = "({0}.{1})".format(neighbor_color, bond_type)
                    color_elements[component] += 1
            for name in sorted(color_elements):
            	atom_color += "({0}_{1})".format(color_elements[name], name) if color_elements[name] > 1 else "({0})".format(name)
            atom_color_with_neighbors[index] = atom_color
        return atom_color_with_neighbors

    def first_round_color(self, isotope_resolved=False, bond_stereo=False, charge=False, R=False, atom_stereo=False, resonance=True, depth=5000):
    	"""
        To do the first round of coloring this compond. We add neighbors' information layer by layer to the atom color identifier until 
        it has a unique identifier or all the atoms in the compound have been used for naming. (based on the readth first search algorithm)
    
        :param isotope_resolved: If true, add isotope information when constructing colors.
        :type isotope_resolved: :py:obj:`True` or :py:obj:`False`.
        :param bond_stereo:  If true, add stereo information to bonds when constructing colors.
        :type bond_stereo: :py:obj:`True` or :py:obj:`False`.
        :param charge: If true, add charge information when constructing colors.
        :type charge: :py:obj:`True` or :py:obj:`False`.
        :param R: If true, add R group information when constructing colors.
        :type R: :py:obj:`True` or :py:obj:`False`.
        :param atom_stereo: If true, add atom stereochemistry information when constructing colors.
        :type atom_stereo: :py:obj:`True` or :py:obj:`False`.
        :param resonance: If true, detect resonant compound pairs without disginguish between double bonds and single bonds.
        :type resonance: :py:obj:`True` or :py:obj:`False`.
        :param int depth: the max depth of coloring.

        :return:
    	"""
        excluded_index = self.metal_index + self.H_index
        excluded_index += self.R_groups if R else []
        atoms_to_color = [i for i in range(len(self.atoms)) if i not in excluded_index]

        self.reset_color()
        self.generate_atom_zero_layer_color(isotope_resolved=isotope_resolved, charge=charge, atom_stereo=atom_stereo)
        atom_color_with_neighbors = self.generate_atom_color_with_neighbors(atoms_to_color, excluded=excluded_index,
        	zero_layer=True, resonance=resonance, bond_stereo=bond_stereo)

        atom_neighbors = collections.defaultdict(list)
        visited = collections.defaultdict(set)
        
        for index in atoms_to_color:
            atom_neighbors[index].append(index)
            visited[index].add(index)
            
        i = 0
        if depth == 5000:
            depth = len(atoms_to_color)

        while i < depth and atoms_to_color:
            current_layer_color_groups = collections.defaultdict(list)
            for atom_index in atoms_to_color:
                color_elements = collections.Counter()
                for neighbor_index in atom_neighbors[atom_index]:
                    color_elements[atom_color_with_neighbors[neighbor_index]] += 1
                added = ""
                for name in sorted(color_elements):
            		added += "({0}_{1})".format(color_elements[name], name) if color_elements[name] > 1 else "({0})".format(name)
                atom.color += added
                atom.color_layers[i+1] = atom.color
                atom_neighbors[atom.atom_number] = self.get_next_layer_neighbors(atom_neighbors[atom_index], visited[atom_index])
                current_layer_color_groups[atom.color].append(atom_index)

            if i > 3:
            	# avoid early stop
                atom_to_color_update = []
                for name in current_layer_color_groups.keys():
                    if len(current_layer_color_groups[name]) > 1:
                        atom_to_color_update.extend(current_layer_color_groups[name])
                atoms_to_color = atom_to_color_update
            i += 1
    
    def invalid_symmetric_atoms(self, bond_stereo=False, resonance=True, R=False):
    	"""
        Check if the atoms with the same color identifier are symmetric.

        :return: the dictionary of atom color identifier with the corresponding invalid symmetric atom index.
    	"""

        excluded_index = self.metal_index + self.H_index
        excluded_index += self.R_groups if R else []
        atoms_to_color = [i for i in range(len(self.atoms)) if i not in excluded_index]

        atom_color_with_neighbors = self.generate_atom_color_with_neighbors(atoms_to_color, excluded=excluded_index, 
        	zero_layer=False, resonance=resonance, bond_stereo=bond_stereo)

        not_valid = []
        color_groups = self.color_groups(excluded=excluded_index)
        for name in color_groups.keys():
            if len(color_groups[name]) > 1:
                atom_index_with_same_color = color_groups[name]
                visited = collections.defaultdict(set)
                for atom_index in atom_index_with_same_color:
                	visited[atom_index].add(atom_index)
                current_layer = {atom_index: [atom_index] for atom_index in atom_index_with_same_color}
                count = 0
                flag = True
                while flag:
                    count += 1
                    target_index = atom_index_with_same_color[0]
                    target_color_list = [atom_color_with_neighbors[atom_index] for atom_index in current_layer[target_index]]
                    target_color_list.sort()
                    target_color = ''.join(target_color_list)

                    "check if neighors of this layer are the same among all the atoms with the same color identifier."
                    for compared_index in current_layer.keys():
                        if compared_index != target_index:
                            compared_color_list = [atom_color_with_neighbors[atom_index] for atom_index in current_layer[compared_index]]
                            compared_color_list.sort()
                            compared_atom_color = ''.join(compared_color_list)

                            if target_color != compared_atom_color:
                                not_valid.append(atom_index_with_same_color)
                                flag = False
                                break

                    "If the share the same neighbors, check the next layer."
                    if flag:
                        for atom_index in current_layer.keys():
                            next_layer_neighbors = self.get_next_layer_neighbors(current_layer[atom_index], visited[atom_index])
                            current_layer[atom_index] = next_layer_neighbors

                        # pruning check, to see if the have the same number of next layer's neighbors.
                        target_length = len(current_layer[target_index])
                        for atom_index in current_layer.keys():
                            if len(current_layer[atom_index]) != target_length:
                                not_valid.append(atom_index_with_same_color)
                                flag = False
                                break

                        # we have checked all the neighbors.
                        if target_length == 0:
                            flag = False
        return not_valid

    def curate_invalid_symmetric_atoms(self, bond_stereo=False, resonance=True):
    	"""
        Curate the atom color identifier of invalid symmetric atom.
        We recolor those invalid atoms with the full color identifiers of its neighbors layer by layer until where the difference can be captured.

        :return:
    	"""
        not_valid = self.invalid_symmetric_atoms()
        while not_valid:
	        atom_color_with_neighbors = self.generate_atom_color_with_neighbors(zero_layer=False, resonance=resonance, bond_stereo=bond_stereo)
	        for invalid_symmetric_atom_index in not_valid:
	            invalid_symmetric_atom_index
	            for atom_index in invalid_symmetric_atom_index:
	                visited = {atom_index}
	                current_layer_neighbors = [atom_index]
	                atom_color = ""
	                while current_layer_neighbors:
	                	color_elements = collections.Counter()
                		for neighbor_index in current_layer_neighbors:
                    		color_elements[atom_color_with_neighbors[neighbor_index]] += 1
	                    added = ""
                		for name in sorted(color_elements):
            				added += "({0}_{1})".format(color_elements[name], name) if color_elements[name] > 1 else "({0})".format(name)
	                    atom_color += added
	                    current_layer_neighbors = self.get_next_layer_neighbors(current_layer_neighbors, visited)
	                self.atoms[atom_index].color = atom_color
        	not_valid = self.invalid_symmetric_atoms()

    def color_metal(self, charge=False, isotope_resolved=False, bond_stereo=False, resonance=True):
    	"""
        To color the metals in the compound. Here we just incorprate information of directly connected atoms.
        :param list atom_index: the list of atom index that needs coloring.
        :param bond_stereo:  If true, add stereo information to bonds when constructing colors.
        :type bond_stereo: :py:obj:`True` or :py:obj:`False`.
        :param charge: If true, add charge information when constructing colors.
        :type charge: :py:obj:`True` or :py:obj:`False`.
        :param isotope_resolved: If true, add isotope information when constructing colors.
        :type isotope_resolved: :py:obj:`True` or :py:obj:`False`.
        :param resonance: If true, detect resonant compound pairs without disginguish between double bonds and single bonds.
        :type resonance: :py:obj:`True` or :py:obj:`False`.

        :return:
    	"""
    	atom_color_with_neighbors = self.generate_atom_color_with_neighbors(self.metal_index, excluded=self.H_index, 
    		zero_layer=False, resonance=False, bond_stereo=False)
    	for atom_index in self.metal_index:
    		self.atoms[atom_index] = atom_color_with_neighbors[atom_index]
      
    def color_H(self, charge=False, isotope_resolved=False, bond_stereo=False, resonance=True):  
    	"""To color the metals in the compound. Here we just incorprate information of directly connected atoms.
        :param list atom_index: the list of atom index that needs coloring.
        :param bond_stereo:  If true, add stereo information to bonds when constructing colors.
        :type bond_stereo: :py:obj:`True` or :py:obj:`False`.
        :param charge: If true, add charge information when constructing colors.
        :type charge: :py:obj:`True` or :py:obj:`False`.
        :param isotope_resolved: If true, add isotope information when constructing colors.
        :type isotope_resolved: :py:obj:`True` or :py:obj:`False`.
        :param resonance: If true, detect resonant compound pairs without disginguish between double bonds and single bonds.
        :type resonance: :py:obj:`True` or :py:obj:`False`.

        :return:
    	"""
    	atom_color_with_neighbors = self.generate_atom_color_with_neighbors(self.H_index, zero_layer=False, 
    		resonance=False, bond_stereo=False)
    	for atom_index in self.H_index:
    		self.atoms[atom_index] = atom_color_with_neighbors[atom_index]

    def metal_color_identifier(self, details=True):
    	"""
        To generate the metal color string representation.
        :param details: if true, use full metal color when constructing identifier.
        :type details: :py:obj:`True` or :py:obj:`False`.

        :return: the metal color string representation.
    	"""
        if details:  
            color_counter = collections.Counter([self.atoms[index].color for index in self.metal_index])
        else:
            color_counter = collections.Counter([self.atoms[index].color_0 for index in self.metal_index])
        return "".join(["({0})({1})".format(color_counter[key], key) for key in sorted(color_counter)])
    
    def H_color_identifier(self, details=True):
    	"""
        To generate the H color string representation.
        :param details: if true, use full H color when constructing identifier.
        :type details: :py:obj:`True` or :py:obj:`False`.

        :return: the H color string representation.
    	"""
        if details:
            color_counter = collections.Counter([self.atoms[index].color for index in self.H_index])
        else:
            color_counter = collections.Counter([self.atoms[index].color_0 for index in self.H_index])
        return "".join(["({0})({1})".format(color_counter[key], key) for key in sorted(color_counter)])

    @property
    def backbone_color_identifier(self):
    	"""
        To generate the backbone color identifier for this compound. Exclude Hs and metals.

        :return: the color string identifier for this compound.
    	"""
    	color_groups = self.color_groups(excluded=self.metal_index + self.H_index)
    	return "".join(["({0})({1})".format(len(color_groups[key]), key) for key in sorted(color_groups)])
