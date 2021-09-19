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
    	valence, h0designator, atom_mappings, inversion_retention_flag, atom_number, exact_change_flag, default_symbol, kat=""):
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
	:param str atom_mappings: atom_mappings.
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
        self.atom_mappings = atom_mappings.strip()
        self.inversion_retention_flag = inversion_retention_flag.strip()
        self.atom_number = atom_number
        self.exact_change_flag = exact_change_flag.strip()
        self.default_symbol = default_symbol.strip()
        self.neighbors = []
        self.color_0 = ""
        self.color = ""
        self.in_cycle = False
        self.bond_counts = 0
        self.group_id = 0
        self.double_bond_counts = 0
        self.color_tuple = tuple()
        self.color_layers = {}
        self.distance_to_R = 0
        self.kat = kat
        self.atom_mappings = []

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
    
    def is_R(self):
    	"""
	To determine if the atom is an R group.

	:return: boolean.
    	"""
        if "A" in self.default_symbol or "R" in self.default_symbol or "*" in self.default_symbol:
            if self.default_symbol not in notR:
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
        self.bond_type = bond_type
        return self.bond_type

    def update_stereochemistry(self, stereo):
    	"""
	To update the stereochemistry of the bond.
	:param str stereo: the updated atom stereochemistry.

	:return: the updated stereochemistry.
    	"""
        self.bond_stereo = stereo
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
        self.color_groups = collections.defaultdict(list)
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
                counter[atom.atom_symbol] += 1
        
        return "".join([ char + str(counter[char]) for char in sorted(counter.keys()) ])
    
    @property
    def R_groups(self):
    	"""
	To get all the R groups in the compound.

	:return: the list of index of all the R groups.
    	"""
        rs = []
        for index, atom in enumerate(self.atoms):
            if atom.is_R():
                rs.append(index)
        return rs

    def contains_R_groups(self):
    	"""
	To check if the compound contains R group(s).

	:return: boolean.
    	"""
        return self.R_groups != []

    @property
    def isolated_metal_index(self):
        """
	To get the indexes of metals without connections to other atoms in the compound.

	:return: the list of index of isolated metals.
        """
        return [index for index, atom in enumerate(self.atoms) if atom.atom_symbol in metalSymbols]
    
    @property
    def HIndex(self):
    	"""
	To get the index of H(s) in the compound. 

	:return: the list of index of Hs.
    	"""
        return [index for index, atom in enumerate(self.atoms) if atom.atom_symbol == "H"]

    @property
    def heavy_atoms(self):
       """
   	To get all the heavy atoms in the compound.

   	:return: the list of heavy atoms in the compound.
       """
        return [atom for atom in self.atoms if atom.atom_symbol != "H"]
    
    @property
    def index_of_heavy_atoms(self):
    	"""
	To map the atom index to heavy atoms.

	:return: the dictionay of atom index to atom entity.
    	"""
        return { self.heavy_atoms[i].atom_number: i for i in range(len(self.heavy_atoms)) }

    def metals(self):
    	"""
	To get all the metals in the compound.

	:return: the dictionary of metal and the corresponding list of index.
    	"""
        metals = collections.defaultdict(list)
        for index, atom in enumerate(self.atoms):
            if atom.default_symbol in metalSymbols:
                metals[atom.default_symbol].append(index)
        return metals

    def update_double_bond_stereo(self, bonds):
	"""
	To update the stereochemistry of the double bonds.
	:param list bonds: the list of bonds whose bond stereo needs update. 

	:return: 
	"""
	for bond in bonds:
	    first_atom = bond.first_atom_number
	    second_atom = bond.second_atom_number
	    if self.bond_lookup[(first_atom, second_atom)].bond_type == "2":
	        self.bond_lookup[(first_atom, second_atom)].bond_stereo = str(bond.bond_stereo)

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
                bond.bond_type = "4"
        bond_out_of_cycle = self.extract_double_bond_connecting_cycle(atom_in_cycle)
        for bond in bond_out_of_cycle:
            bond.bond_type = "1"


    def extract_double_bond_connecting_cycle(self, atom_in_cycle):
	"""
	Extract the double bonds connecting to the aromatic cycles.
	:param list atom_in_cycle: the list of indexes of atoms in the aromatic cycles.

	:return: the list of double bonds connecting to the aromatic cycles.
	"""
        double_bond_connecting_cycle = []
        for atom_index in atom_in_cycle:
            if self.atoms[atom_index].atom_symbol == "C":
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
                        if self.atoms[neighbor].atom_symbol != "H":
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

    def color_compound(self, R=True,stereo_resolved=True, atom_stereo=True, strict=True):

        self.first_round_color(R=R, stereo_resolved=stereo_resolved, atom_stereo=atom_stereo, strict=strict)
        self.curate_symmetry(stereo_resolved=stereo_resolved, strict=strict)
        self.color_metals(stereo_resolved=stereo_resolved, strict=strict)
    
    def update_color_tuple(self, resonance=False):
    	"""
        To update the color tuple of the atoms in the compound. This color tuple includes information of its neighboring atoms and bonds.
	:param boolean resonance: to find the resonant structure. This will ignore the difference between single and double bonds.

        :return:
    	"""
        for atom in self.heavy_atoms:
            elements = collections.Counter()
            bond_types = collections.Counter()
            for neighbor_index in atom.neighbors:
                neighbor = self.atoms[neighbor_index]
                if neighbor.atom_symbol != "H":
                    elements[neighbor.atom_symbol] += 1
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
	param the_other: the mappings compound entity.
	:type the_other: :class:`~MDH.compound.Compound`
	:param boolean resonance: to find the resonant structure. This will ignore the difference between single and double bonds.
	:param boolean R: to take account of positon of R groups.
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

    def match_loosely(self, the_other, R=False):
    	"""

    	"""
    	mappings = self.find_mappings(self, the_other, resonance=True, R=False)
        one_to_one_mappings = self.generate_one_to_one_mappings(the_other, mappings)

        validMatch = []
        for cur_mappings in one_to_one_mappings:
            flag = True
            reversed_mappings = { cur_mappings[key] : key for key in cur_mappings }
            for from_index in sorted(cur_mappings.keys()):
                to_index = cur_mappings[from_index]
                if 1 in self.atoms[from_index].color_layers and 1 in the_other.atoms[to_index].color_layers and 
                self.atoms[from_index].color_layers[1] == the_other.atoms[to_index].color_layers[1]:
                    continue
                
                three = {from_index}
                from_index_neighbor = self.find_double_bond_linked_atom(from_index)
                to_index_neighbor = the_other.find_double_bond_linked_atom(to_index)
                
                if from_index_neighbor == -1 and to_index_neighbor == -1:
                    flag = False
                    break

                elif from_index_neighbor == -1:
                    reversed_from_index_neighbor = reversed_mappings[to_index_neighbor]
                    three.add(reversed_from_index_neighbor)
                    fNextNeighbor = self.find_double_bond_linked_atom(reversed_from_index_neighbor)
                    if fNextNeighbor != -1:
                        three.add(fNextNeighbor)

                elif to_index_neighbor == -1 and from_index_neighbor in cur_mappings:
                    three.add(from_index_neighbor)
                    reverseto_index_neighbor = cur_mappings[from_index_neighbor]
                    nNextNeighbor = the_other.find_double_bond_linked_atom(reverseto_index_neighbor)
                    if nNextNeighbor != -1:
                        three.add(reversed_mappings[nNextNeighbor])

                elif from_index_neighbor != -1 and to_index_neighbor != -1:
                    three.add(from_index_neighbor)
                    three.add(reversed_mappings[to_index_neighbor])
                
                if len(three) < 3:
                    flag = False
                    break
                
                k = [self.atoms[i].atom_symbol for i in three].count("C")
                if k == 3:
                    flag = False
                    break
            if flag:
                validMatch.append(cur_mappings)
        return []
    
    def find_double_bond_linked_atom(self, i):
    	"""
        Find the index of the atom that is doubly linked to the target atom[i].
        param int i: the index of the target atom.

        :return: the index of the doubly linked atom.
    	"""
        for neighbor_index in self.atoms[i].neighbors:
            if neighbor_index in self.index_of_heavy_atoms and self.bond_lookup[(i, neighbor_index)].bond_type == "2":
                return neighbor_index
        return -1

    def detect_linear_circular_interchange(self, the_other):
    	"""
	Detect if the two compounds 
	:param the_other: the mappings compound entity.
	:type the_other: :class:`~MDH.compound.Compound`

	:return: boolean.
    	"""
        sAtoms = [atom.color_layers[1] if 1 in atom.color_layers else atom.atom_symbol for atom in self.heavy_atoms]
        aAtoms = [atom.color_layers[1] if 1 in atom.color_layers else atom.atom_symbol for atom in self.heavy_atoms]
        atom_1ntersect = list((collections.Counter(sAtoms) & collections.Counter(aAtoms)).elements())
        #print("atom name cnt: ", len(sAtoms), len(aAtoms), len(atom_1ntersect))
        if 2 <= len(sAtoms) - len(atom_1ntersect) <= 3:
            sinCycle = [ atom.in_cycle for atom in self.heavy_atoms].count(True)
            ainCycle = [ atom.in_cycle for atom in the_other.heavy_atoms].count(True)
            #print("incycle number: ", sinCycle, ainCycle)
            if sinCycle != ainCycle:
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
		\___/		\___/
		 ___         ___
		/	\       /   \ 
	   L 	 H     L     L
	   	trans         cis
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

        if (first_light_side and first_light_side.atom_number == first_heavy_side.atom_number) or 
        (second_light_side and second_light_side.atom_number == second_heavy_side.atom_number):
        	return 0
       
        if vertical:
            if (first_heavy_side.x - first_atom.x) * (second_heavy_side.x - second_atom.x) > 0:
                return 1
            else:
                return -1
        else:
            if (self.calculate_y_coordinate(slope, b, first_heavy_side) - first_heavy_side.y) * (self.calculate_y_coordinate(slope, b, second_heavy_side) - second_heavy_side.y) > 0:
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

    def rank_sides(self, neighbors, atom_forming_double_bond):
    	"""
    To determine the rank of the two branches connecting the atom forming the double bond. This is based on comparison of the atomic weights of the two branches.
    Breadth first algorithm.
    :param list neighbors: the list of atom indexed of the atoms connecting the atom forming the double bond.
    :param atom_forming_double_bond: the bond entity.
	:type atom_forming_double_bond: :class:`~MDH.compound.Atom`
    	"""
    	if len(neighbors) < 2:
    		return self.atoms[neighbors[0]], None

    	one_atom, the_other_atom = self.atoms[neighbors[0]], self.atoms[neighbors[1]]

        if atomic_weights[one_atom.default_symbol] > atomic_weights[the_other_atom.default_symbol]:
            return [one_atom, the_other_atom]

        elif atomic_weights[one_atom.default_symbol] < atomic_weights[the_other_atom.default_symbol]:
            return [the_other_atom, one_atom]

        else:
            one_neighbors = [(one_atom, atom_forming_double_bond.atom_number)]
            one_visited = set([atom_forming_double_bond.atom_number])
            the_other_neighbors = [(the_other_atom, atom_forming_double_bond.atom_number)]
            the_other_visited = set([patom_forming_double_bond.atom_number])
            one_neighbor_atomic_weight_list = self.collect_atomic_weights_of_neighbors(one_neighbors)
            the_other_neighbor_atomic_weight_list = self.collect_atomic_weights_of_neighbors(the_other_neighbors)

            while one_neighbor_atomic_weight_list == the_other_neighbor_atomic_weight_list and one_neighbor_atomic_weight_list:
                one_neighbors_updated = []
                for atom, pre_index in one_neighbors:
                    one_visited.add(atom.atom_number)
                    for neighbor in atom.neighbors:
                        if neighbor != pre_index and neighbor not in one_visited:
                            one_neighbors_updated.append((self.atoms[neighbor], atom.atom_number))
                one_neighbors = one_neighbors_updated
                one_neighbor_atomic_weight_list = self.collect_atomic_weights_of_neighbors(one_neighbors)

                the_other_neighbors_updated = []
                for atom, pre_index in the_other_neighbors:
                    the_other_visited.add(atom.atom_number)
                    for neighbor in atom.neighbors:
                        if neighbor!= pre_index and neighbor not in the_other_visited:
                            the_other_neighbors_updated.append((self.atoms[neighbor], atom.atom_number))
                the_other_neighbors = the_other_neighbors_updated
                the_other_neighbor_atomic_weight_list = self.collect_atomic_weights_of_neighbors(the_other_neighbors)


            if tuple(one_neighbor_atomic_weight_list) > tuple(the_other_neighbor_atomic_weight_list):
                return [one_atom, the_other_atom]
            elif tuple(one_neighbor_atomic_weight_list) < tuple(the_other_neighbor_atomic_weight_list):
                return [the_other_atom, one_atom]
            else:
                return [one_atom, one_atom]

    def collect_atomic_weights_of_neighbors(self, neighbors):
    	"""
	To collect the atomic weights of the next layer's neighbors.
	:param list neighbors: the list of this layer's atom with its parent atom. The parent is used to avoid recalculating the former layer.

	:return: the list of atomic weights of next layer's neighbors
    	"""
        updated_neighbors = []
        for atom, pre in neighbors:
            updated_neighbors.extend([atomic_weights[self.atoms[neighbor].default_symbol] for neighbor in atom.neighbors if neighbor != pre])
        updated_neighbors.sort(reverse=True)
        return updated_neighbors

    def isolatedCpd(self):

        for atom in self.atoms:
            if atom.bond_counts == 0:
                return True
        return False

    

    def abnormalBond(self):

        abnormalAtoms = collections.defaultdict(list)
        for index, atom in enumerate(self.atoms):
            if atom.atom_symbol in atomBonds:
                bond_counts = atom.bond_counts
                bond_counts -= atom.charge
                if type(atomBonds[atom.atom_symbol]) is list:
                    if bond_counts not in atomBonds[atom.atom_symbol]:
                        abnormalAtoms[atom.atom_symbol].append(index)
                else:
                    if bond_counts > atomBonds[atom.atom_symbol]:
                        abnormalAtoms[atom.atom_symbol].append(index)
        return abnormalAtoms

    def curateInvalidN(self):

        abnormalAtoms = self.abnormalBond()
        if "N" in abnormalAtoms:
            for atom_index in abnormalAtoms["N"]:
                atom = self.atoms[atom_index]
                atom.charge += 1
                for neighbor in atom.neighbors:
                    if self.atoms[neighbor].atom_symbol == "O" and self.bond_lookup[(atom.atom_number, neighbor)].bond_type == "2":
                        self.atoms[neighbor].charge -= 1
                        self.bond_lookup[(atom.atom_number, neighbor)].bond_type = "1"

    def clean_colors(self):

        for atom in self.atoms:
            atom.color_0 = ""
            atom.color = ""
            atom.color_layers = {}
        self.color_groups = collections.defaultdict(list)
    


    def first_round_color(self, isotope_resolved=False, stereo_resolved=False, charge=False, R=False, atom_stereo=False, strict=True, depth=5000):
        
        atomNeighbors = collections.defaultdict(list)
        color_groups = collections.defaultdict(list)
        excluded = self.isolated_metal_index + self.HIndex
        self.clean_colors()

        for index, atom in enumerate(self.atoms):
            if index not in excluded:
                if not atom.is_R():
                    atom.color_0 += atom.atom_symbol
                else:
                    if R:
                        atom.color_0 += "R"
                    else:
                        continue
                if atom_stereo:
                    atom.color_0 += atom.atom_stereo_parity
                if charge:
                    atom.color_0 += str(atom.charge)
                if isotope_resolved:
                    atom.color_0 += atom.mass_difference
                atom.color = atom.color_0
                atomNeighbors[index].append(atom)
                color_groups[atom.color].append(atom)

        
        atomNeighborRepresentation = collections.defaultdict(str)
        numOfValidAtom = 0

        notComplete = []
        for index, atom in enumerate(self.atoms):
            if atom.color_0:
                color_elements = []
                numOfValidAtom += 1
                for neighbor_index in atom.neighbors:
                    if self.atoms[neighbor_index].color_0:
                        connecting_bond = self.bond_lookup[(atom.atom_number, neighbor_index)]
                        bond_type = connecting_bond.bond_type
                        if not strict and bond_type == "2":
                            bond_type = "1"
                        if stereo_resolved:
                            color_elements.append((self.atoms[neighbor_index].color_0, bond_type + connecting_bond.bond_stereo))
                        else:
                            color_elements.append((self.atoms[neighbor_index].color_0, bond_type))
                notComplete.append(atom)
                atomNeighborRepresentation[index] = atom.color_0 + "".join(["(" + a + "," + b + ")" for (a, b) in sorted(color_elements)])
        
        visited = collections.defaultdict(set)
        i = 0
        if depth != 5000:
            numOfValidAtom = depth

        while i < numOfValidAtom and notComplete:
        
            inner_color_groups = collections.defaultdict(list)
            for atom in notComplete:
                newNeighbors = []
                color_list = []
                for neighbor in atomNeighbors[atom.atom_number]:
                    color_list.append(atomNeighborRepresentation[neighbor.atom_number])
                    if neighbor.atom_number not in visited[atom.atom_number]:
                        visited[atom.atom_number].add(neighbor.atom_number)
                        newNeighbors.extend([self.atoms[newNeighbor_index] for newNeighbor_index in neighbor.neighbors if self.atoms[newNeighbor_index].color_0])
                if color_list:
                    color_list.sort()
                    atom.color += "(" + "".join(color_list) + ")"
                    atom.color_layers[i+1] = atom.color
                atomNeighbors[atom.atom_number] = newNeighbors
                inner_color_groups[atom.color].append(atom)
            if i > 3:
                newNotComplete = []
                for name in inner_color_groups.keys():
                    if len(inner_color_groups[name]) > 1:
                        newNotComplete.extend(inner_color_groups[name])
                notComplete = newNotComplete
            i += 1
        
        for atom in self.atoms:
            if atom.color_0:
                self.color_groups[atom.color].append(atom.atom_number)
    
    
    def check_symmetry(self, stereo_resolved=False, strict=True):

        atomNeighborRepresentation = collections.defaultdict(str)
        for atom in self.atoms:
            if atom.color:
                color_elements = []
                for neighbor_index in atom.neighbors:
                    if self.atoms[neighbor_index].color:
                        connecting_bond = self.bond_lookup[(atom.atom_number, neighbor_index)]
                        bond_type = connecting_bond.bond_type
                        if not strict and bond_type == "2":
                            bond_type = "1"
                        if stereo_resolved:
                            color_elements.append((self.atoms[neighbor_index].color, bond_type + connecting_bond.bond_stereo))
                        else:
                            color_elements.append((self.atoms[neighbor_index].color, bond_type))
                atomNeighborRepresentation[atom.atom_number] = atom.color + "".join(["(" + a + "," + b + ")" for (a, b) in sorted(color_elements)])
        
        notValid = {}
        for name in self.color_groups.keys():
            if len(self.color_groups[name]) > 1:
                sameColors = self.color_groups[name]
                
                visited = collections.defaultdict(set)
                currentLayer = {atom_index: [atom_index] for atom_index in sameColors}
                count = 0
                flag = True
                while flag:
                    count += 1
                    targetatom_index = sameColors[0]
                    targetColorList = [atomNeighborRepresentation[atom_index] for atom_index in currentLayer[targetatom_index]]
                    targetColorList.sort()
                    targetColoring = ''.join(targetColorList)

                    for compareatom_index in currentLayer.keys():
                        if compareatom_index != targetatom_index:
                            compareColorList = [atomNeighborRepresentation[atom_index] for atom_index in currentLayer[compareatom_index]]
                            compareColorList.sort()
                            compareColoring = ''.join(compareColorList)

                            if targetColoring != compareColoring:
                                notValid[name] = (sameColors, count)
                                flag = False
                                break
                    if flag:
                        for atom_index in currentLayer.keys():
                            nextNeighbors = []
                            for neighbor_index in currentLayer[atom_index]:
                                if neighbor_index not in visited[atom_index]:
                                    visited[atom_index].add(neighbor_index)
                                    nextNeighbors.extend([newIndex for newIndex in self.atoms[neighbor_index].neighbors if self.atoms[newIndex].color])
                            currentLayer[atom_index] = nextNeighbors
                        targetLength = len(currentLayer[targetatom_index])
                        for atom_index in currentLayer.keys():
                            if len(currentLayer[atom_index]) != targetLength:
                                notValid[name] = (sameColors, count)
                                flag = False
                                break
                        if targetLength == 0:
                            flag = False
        return notValid

    def curate_symmetry(self, stereo_resolved=False, strict=True):

        notValid = self.check_symmetry()
        if not notValid:
            return
        
        atomNeighborRepresentation = collections.defaultdict(str)

        for atom in self.atoms:
            if atom.color:
                color_elements = []
                for neighbor_index in atom.neighbors:
                    if self.atoms[neighbor_index].color:
                        connecting_bond = self.bond_lookup[(atom.atom_number, neighbor_index)]
                        bond_type = connecting_bond.bond_type
                        if not strict and bond_type == "2":
                            bond_type = "1"
                        if stereo_resolved:
                            color_elements.append((self.atoms[neighbor_index].color_0, bond_type + connecting_bond.bond_stereo))
                        else:
                            color_elements.append((self.atoms[neighbor_index].color_0, bond_type))
                atomNeighborRepresentation[atom.atom_number] = atom.color + "".join(["(" + a + "," + b + ")" for (a, b) in sorted(color_elements)])

        
        for name in notValid:
            del self.color_groups[name]

            notSameIndex, layer = notValid[name]
            
            for targetatom_index in notSameIndex:
                visited = set()
                currentLayer = [targetatom_index]
                targetColor = ""
                count = layer + 1
                while count  and currentLayer:
                    curColorList = [atomNeighborRepresentation[atom_index] for atom_index in currentLayer]
                    curColorList.sort()
                    targetColor += "".join(curColorList)
                    nextNeighbors = []
                    for neighbor_index in currentLayer:
                        if neighbor_index not in visited:
                            visited.add(neighbor_index)
                            nextNeighbors.extend([newIndex for newIndex in self.atoms[neighbor_index].neighbors if self.atoms[newIndex].color])
                    currentLayer = nextNeighbors
                    count -= 1
                self.color_groups[targetColor].append(targetatom_index)
                self.atoms[targetatom_index].color = targetColor

        if self.check_symmetry():
            print(self.compund_name)
      
    def color_metals(self, charge=False, isotope_resolved=False, stereo_resolved=False, strict=True):
        
        for atom_index in self.isolated_metal_index:
            atom = self.atoms[atom_index]
            colorelements = collections.defaultdict(int)
            for neighbor_index in atom.neighbors:
                if self.atoms[neighbor_index].color:
                    connecting_bond = self.bond_lookup[(atom.atom_number, neighbor_index)]
                    bond_type = connecting_bond.bond_type
                    if not strict and bond_type == "2":
                        bond_type = "1"
                    if stereo_resolved:
                        component = "({0}, {1})".format(self.atoms[neighbor_index].color, bond_type+connecting_bond.bond_stereo)
                    else:
                        component = "({0}, {1})".format(self.atoms[neighbor_index].color, bond_type)
                    colorelements[component] += 1
            atom.color_0 = atom.atom_symbol
            if charge:
                atom.color_0 += str(atom.charge)
            if isotope_resolved:
                atom.color_0 += atom.mass_difference
            atom.color = atom.color_0
            for name in sorted(colorelements.keys()):
                atom.color += "({0})({1})".format(colorelements[name], name)

    def color_H(self, stereo_resolved=False, charge=False, isotope_resolved=False):

        for index in self.HIndex:
            atom = self.atoms[index]
            colorelements = []
            for neighbor_index in atom.neighbors:
                connecting_bond = self.bond_lookup[(atom.atom_number, neighbor_index)]
                if stereo_resolved:
                    colorelements.append((self.atoms[neighbor_index].color, connecting_bond.bond_type + connecting_bond.bond_stereo))
                else:
                    colorelements.append((self.atoms[neighbor_index].color, connecting_bond.bond_type))
            colorelements.sort()
            atom.color_0 = atom.atom_symbol
            if charge:
                atom.color_0 += str(atom.charge)
            if isotope_resolved:
                atom.color_0 += atom.mass_difference
            atom.color = atom.color_0 + "".join(["(" + a + "," + b + ")" for (a, b) in sorted(colorelements)])


    
   
    def metalString(self, noBond=True):
        if noBond:  
            stringCounter = collections.Counter([self.atoms[index].color_0 for index in self.isolated_metal_index])
        else:
            stringCounter = collections.Counter([self.atoms[index].color for index in self.isolated_metal_index])
        stringRepresent = ""
        #print(stringCounter)
        keys = sorted(stringCounter)
        for key in keys:
            stringRepresent += "({0})({1})".format(stringCounter[key], key)
        return stringRepresent
    
    def HString(self, details=False):
        if details:
            stringCounter = collections.Counter([self.atoms[index].color for index in self.HIndex])
        else:
            stringCounter = collections.Counter([self.atoms[index].color_0 for index in self.HIndex])
        stringRepresent = ""
        keys = sorted(stringCounter)
        for key in keys:
            stringRepresent += "({0})({1})".format(stringCounter[key], key)
        return stringRepresent

    @property
    def colorMetalNoBond(self):
        return self.backBoneString + self.metalString()

    @property
    def colorMetal(self):
        return self.backBoneString + self.metalString(noBond=False)
    
    @property
    def backBoneString(self):
        stringRepresent = ""
        keys = sorted(self.color_groups)
        for key in keys:
            stringRepresent += "({0})({1})".format(len(self.color_groups[key]), key)
        return stringRepresent

    

   

    
    

