#!/usr/bin/python3


"""
MDH.compound
~~~~~~~~~~~~

This module provides the :class:`~MDH.compound.Atom` class, the :class:`MDH.compound.Bond` class,
and the :class:`~MDH.compound.Compound` class to construct compound entity.

"""

import collections
import numpy
import heapq
from pathlib import Path
import ctfile
from . import BASS
from .supplement import not_r_groups
from .supplement import standard_bond_counts
from .supplement import atomic_weights
from .supplement import metal_symbols
from .supplement import index_to_charge
# from .supplement import charge_to_index

class Atom:

    """ Atom class describes the :class:`~MDH.compound.Atom` entity in the compound. """

    def __init__(self, atom_symbol, atom_number, x="0", y="0", z="0", mass_difference="0", charge="0", atom_stereo_parity="0",
                 hydrogen_count="0", stereo_care_box="0", valence="0", h0designator="0", atom_atom_mapping_number="0",
                 inversion_retention_flag="0", exact_change_flag="0", kat=""):
        """Atom initializer.

        :param atom_symbol: atom_symbol.
        :type atom_symbol: :py:obj:`str`.
        :param atom_number: atom_number.
        :type atom_number: :py:obj:`int`.
        :param x: the atom x coordinate.
        :type x: :py:obj:`str`.
        :param y: the atom y coordinate.
        :type y: :py:obj:`str`.
        :param z: the atom z coordinate.
        :type z: :py:obj:`str`.
        :param mass_difference: difference from mass in periodic table.
        :type mass_difference: :py:obj:`str`.
        :param charge: charge.
        :type charge: :py:obj:`str`.
        :param atom_stereo_parity: atom stereo parity.
        :type atom_stereo_parity: :py:obj:`str`.
        :param hydrogen_count: hydrogen_count.
        :type hydrogen_count: :py:obj:`str`.
        :param stereo_care_box: stereo_care_box.
        :type stereo_care_box: :py:obj:`str`.
        :param valence: valence.
        :type: valence: :py:obj:`str`.
        :param h0designator: h0designator.
        :type h0designator: :py:obj:`str`.
        :param atom_atom_mapping_number: atom_atom_mapping_number.
        :type atom_atom_mapping_number: :py:obj:`str`.
        :param inversion_retention_flag: inversion_retention_flag.
        :type inversion_retention_flag: :py:obj:`str`.
        :param exact_change_flag: exact_change_flag.
        :type exact_change_flag: :py:obj:`str`.
        :param kat: KEGG atom type.
        :type kat: :py:obj:`str`.
    	"""
        self.x = float(x.strip())
        self.y = float(y.strip())
        self.z = float(z.strip())
        self.atom_symbol = atom_symbol.strip()
        self.mass_difference = mass_difference.strip()
        self.charge = index_to_charge[charge.strip()] if charge.strip() in index_to_charge else 0
        self.atom_stereo_parity = atom_stereo_parity.strip() if atom_stereo_parity.strip() != "3" else "0"
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
        self.color_layers = collections.defaultdict()
        self.in_cycle = False
        self.bond_counts = 0
        self.group_id = 0
        self.double_bond_counts = 0
        self.distance_to_r = 0
        self.kat = kat

    def update_symbol(self, symbol):
        """
        To update the atom symbol.

        :param symbol: the updated atom symbol.
        :type symbol: :py:class:`str`.
        :return: the updated atom_symbol.
        :rtype: :py:class:`str`.
    	"""
        self.atom_symbol = symbol
        return self.atom_symbol

    def update_atom_number(self, index):
        """
        To update the atom number.

        :param index: the updated atom number.
        :type index: :py:class:`int`.
        :return: the updated atom number.
        :rtype: :py:class:`int`.
        """
        self.atom_number = index
        return self.atom_number

    @property
    def default_symbol(self):
        """
        To get the atom symbol. If the atom symbol is not specified, return "R".

        :return: the atom symbol.
        :rtype: :py:class:`str`.
        """
        if self.is_r:
            return "R"
        return self.atom_symbol

    def remove_neighbors(self, neighbors):
        """
        To remove neighbors to the atom.

        :param neighbors: the list of neighbors that will be removed from this atom.
        :type neighbors: :py:class:`list`.
        :return: the updated list of neighbors of the atom.
        :rtype: :py:class:`list`.
        """
        for atom_index in neighbors:
            if atom_index in self.neighbors:
                self.neighbors.remove(atom_index)
        return self.neighbors

    def add_neighbors(self, neighbors):
        """
        To add neighbors to the atom.

        :param neighbors: the list of neighbors that will be removed from this atom.
        :type neighbors: :py:class:`list`.
        :return: the updated list of neighbors of the atom.
        :rtype: :py:class:`list`.
        """
        for atom_index in neighbors:
            if atom_index not in self.neighbors:
                self.neighbors.append(atom_index)
        return self.neighbors

    def update_stereochemistry(self, stereo):
        """
        To update the atom stereochemistry.

        :param stereo: the updated atom stereochemistry.
        :type stereo: :py:class:`str`.
        :return: the updated atom stereochemistry.
        :rtype: :py:class:`str`.
    	"""
        self.atom_stereo_parity = stereo
        return self.atom_stereo_parity
    
    def color_atom(self, isotope_resolved=False, charge=False, atom_stereo=False):
        """
        To generate the zero layer atom color.

        :param isotope_resolved: If true, add isotope information when constructing colors.
        :type isotope_resolved: :py:obj:`bool`.
        :param charge: If true, add charge information when constructing colors.
        :type charge: :py:obj:`bool`.
        :param atom_stereo: If true, add atom stereochemistry information when constructing colors.
        :type atom_stereo: :py:obj:`bool`.
	    :return: the zero layer atom color.
	    :rtype: :py:class:`str`.
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
        Reset the atom color.

        :return: None.
        :rtype: :py:obj:`None`.
    	"""
        self.color_0 = ""
        self.color = ""
        self.color_layers = {}
    
    @property
    def is_r(self):
        """
        To determine if the atom symbol is specified.

        :return: bool whether the atom symbol is specified.
        :rtype: :py:obj:`bool`.
        """
        if "A" in self.atom_symbol or "R" in self.atom_symbol or "*" in self.atom_symbol or self.atom_symbol == "X":
            if self.atom_symbol not in not_r_groups:
                return True
        return False

    def update_kat(self, kat):
        """
        To update the atom KEGG atom type.

        :param kat: the KEGG atom type for this atom,
        :type kat: :py:class:`str`.
        :return: the updated KEGG atom type.
        :rtype: :py:class:`str`.
        """
        self.kat = kat
        return self.kat

    def clone(self):
        """
        To clone the atom.

        :return: the cloned atom.
        :rtype: :class:`~MDH.compound.Atom`.
        """
        return Atom(self.atom_symbol, self.atom_number, x=str(self.x), y=str(self.y), z=str(self.z), mass_difference=self.mass_difference,
                    charge=str(self.charge), atom_stereo_parity=self.atom_stereo_parity, hydrogen_count=self.hydrogen_count,
                    stereo_care_box=self.stereo_care_box, valence=self.valence, h0designator=self.h0designator,
                    atom_atom_mapping_number=self.atom_atom_mapping_number, inversion_retention_flag=self.inversion_retention_flag,
                    exact_change_flag=self.exact_change_flag, kat=self.kat)

class Bond:

    """ Bond class describes the :class:`~MDH.compound.Bond` entity in the compound. """

    def __init__(self, first_atom_number, second_atom_number, bond_type, bond_stereo="0", bond_topology="0", reacting_center_status="0"):
        """Bond initializer.

        :param first_atom_number: the index of the first atom forming this bond.
        :type first_atom_number: :py:class:`str`.
        :param second_atom_number: the index of the second atom forming this bond.
        :type second_atom_number: :py:class:`str`.
        :param bond_type: the bond type. (1 = Single, 2 = Double, 3 = Triple, 4 = Aromatic, 5 = Single or
        Double, 6 = Single or Aromatic, 7 = double or Aromatic 8 = Any)
        :type bond_type: :py:class:`str`.
        :param bond_stereo: the bond stereo. (Single bonds: 0 = not stereo, 1 = Up, 4 = Either, 6 = Down;
        Double bonds: determined by x, y, z coordinates)
        :type bond_stereo: :py:class:`str`.
        :param bond_topology: bond topology. (O = Either, 1 = Ring, 2 = Chain)
        :type bond_topology: :py:class:`str`.
        :param reacting_center_status: reacting center status.
        :type reacting_center_status: :py:class:`str`.
    	"""
        self.first_atom_number = int(first_atom_number) - 1
        self.second_atom_number = int(second_atom_number) - 1
        self.bond_type = bond_type.strip()
        self.bond_stereo = bond_stereo.strip() if bond_stereo.strip() != "4" and bond_stereo.strip() != "8" else "0"
        self.bond_topology = bond_topology.strip()
        self.reacting_center_status = reacting_center_status.strip()
 
    def update_bond_type(self, bond_type):
        """
        To update the bond type.

        :param bond_type: the updated bond type.
        :type bond_type: :py:class:`str`.
        :return: the updated bond type.
        :rtype: :py:class:`str`.
    	"""
        self.bond_type = str(bond_type)
        return self.bond_type

    def update_stereochemistry(self, stereo):
        """
        To update the bond stereochemistry.

        :param stereo: the updated bond stereochemistry.
        :type stereo: :py:class:`str`.
        :return: the updated bond stereochemistry.
        :rtype: :py:class:`str`.
    	"""
        self.bond_stereo = str(stereo)
        return self.bond_stereo

    def update_first_atom(self, index):
        """
        To update the first atom number of the bond.

        :param index: the updated first atom number.
        :type index: :py:class:`int`.
        :return: the updated first atom number.
        :rtype: :py:class:`int`.
        """
        self.first_atom_number = index
        return self.first_atom_number

    def update_second_atom(self, index):
        """
        To update the second atom number of the bond.

        :param index: the updated second atom number.
        :type index: :py:class:`int`.
        :return: the updated second atom number.
        :rtype: :py:class:`int`.
        """
        self.second_atom_number = index
        return self.second_atom_number

    def clone(self):
        """
        To clone the bond.

        :return: the cloned bond.
        :rtype: :class:`~MDH.compound.Bond`.
        """
        return Bond(str(self.first_atom_number+1), str(self.second_atom_number+1), self.bond_type, bond_stereo=self.bond_stereo,
                    bond_topology=self.bond_topology, reacting_center_status=self.reacting_center_status)


class Compound:

    """ Compound class describes the :class:`~MDH.compound.Compound` entity. """

    def __init__(self, compound_name, atoms, bonds):
        """Compound initializer.
			
        :param compound_name: the compound name.
        :type compound_name: :py:class:`str`.
        :param atoms: a list of :class:`~MDH.compound.Atom` entities in the compound.
        :type atoms: :py:class:`list`.
        :param bonds: a list of :class:`~MDH.compound.Bond` entities in the compound.
        :type bonds: :py:class:`list`.
    	"""
        self.compound_name = compound_name
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
        self.cycles = self.find_cycles()
        self.calculate_distance_to_r_groups()

    def encode(self):
        """
        To clone the compound.

        :return: the cloned compound.
        :rtype: :class:`~MDH.compound.Compound`.
        """
        return self.compound_name, [atom.clone() for atom in self.atoms], [bond.clone() for bond in self.bonds]

    @property
    def name(self):
        """
        To get the compound name.
        :return: the compound name.
        :rtype: :py:class:`str`.
        """
        return self.compound_name

    @staticmethod
    def create(molfile):
        """
        Create the compound entity based on the molfile representation.

        :param molfile: the filename of the molfile.
        :type molfile: :py:class:`str`.
        :return: the constructed compound entity.
        :rtype: :class:`~MDH.compound.Compound`.
        """
        compound_name = Path(molfile).stem
        try:
            with open(molfile, 'r') as infile:
                ct_object = ctfile.load(infile)
            atoms = [ Atom(atom.atom_symbol, i, atom['x'], atom['y'], atom['z'],  atom['mass_difference'], atom.charge,
                           atom['atom_stereo_parity'], atom['hydrogen_count'], atom['stereo_care_box'], atom['valence'],
                           atom['h0designator'], atom['atom_atom_mapping_number'], atom['inversion_retention_flag'],
                           atom['exact_change_flag']) for i, atom in enumerate(ct_object.atoms) ]
            bonds = [ Bond(bond['first_atom_number'], bond['second_atom_number'], bond['bond_type'], bond['bond_stereo'],
                           bond['bond_topology'], bond['reacting_center_status']) for bond in ct_object.bonds ]
            return Compound(compound_name, atoms, bonds)
        except:
            print(compound_name + "cannot not be converted to ctifle object")
            return None

    @property
    def formula(self):
        """
        To construct the formula of this compound (only consider heavy atoms).

        :return: string formula of the compound.
        :rtype: :py:class:`str`.
    	"""
        counter = collections.Counter()
        for atom in self.atoms:
            if atom.default_symbol == "H":
                continue
            elif atom.is_r:
                counter["R"] += 1
            else:
                counter[atom.default_symbol] += 1
        return "".join([ char + str(counter[char]) for char in sorted(counter.keys()) ])

    @property
    def composition(self):
        """
        To get the atom symbols and bond types in the compound.

        :return: the atom and bond information of the compound
        :rtype: :py:class:`dict`.
        """
        atom_composition = collections.Counter([atom.default_symbol for atom in self.atoms])
        bond_composition = collections.Counter([bond.bond_type for bond in self.bonds])
        atom_composition.update(bond_composition)
        return atom_composition

    @property
    def r_groups(self):
        """
        To get all the R groups in the compound.

        :return: the list of index of all the R groups.
        :rtype: :py:obj:`list`.
    	"""
        return [index for index, atom in enumerate(self.atoms) if atom.is_r]

    def contains_r_groups(self):
        """
        To check if the compound contains R group(s).

        :return: bool whether the compound contains R group.
        :rtype: :py:obj:`bool`.
    	"""
        return self.r_groups != []

    def has_isolated_atoms(self):
        """
        To check if the compound has atoms that have no connections to other atoms.

        :return: bool whether the compound has isolated atoms.
        :rtype: :py:obj:`bool`.
    	"""
        for atom in self.atoms:
            if atom.bond_counts == 0:
                return True
        return False

    @property
    def metal_index(self):
        """
        To get the metals in the compound.

        :return: a list of atom numbers of metals.
        :rtype: :py:class:`list`.
        """
        return [index for index, atom in enumerate(self.atoms) if atom.default_symbol in metal_symbols]
    
    @property
    def h_index(self):
        """
        To get the H(s) in the compound.

        :return: a list of atom numbers of H(s).
        :rtype: :py:class:`list`.
    	"""
        return [index for index, atom in enumerate(self.atoms) if atom.default_symbol == "H"]

    @property
    def heavy_atoms(self):
        """
        To get all the heavy atoms in the compound.

        :return: a list of atom numbers of heavy atoms.
        :rtype: :py:class:`list`.
        """
        return [atom for atom in self.atoms if atom.default_symbol != "H"]

    @property
    def index_of_heavy_atoms(self):
        """
        To map the atom number to index in the heavy atom list.

        :return: the dictionary of atom number to atom index of heavy atoms.
        :rtype: :py:class:`dict`.
    	"""
        return { atom.atom_number: i for i, atom in enumerate(self.heavy_atoms) }

    # def is_symmetric(self, atom_index):
    #
    #     for index in atom_index:
    #         color = self.atoms[index].color
    #         if len(self.color_compound()[color]) > 1:
    #             return False
    #     return True
       
    def color_groups(self, excluded=None):
        """
        To update the compound color groups after coloring.

        :return: the dictionary of atom color with the list of atom number.
        :rtype: :py:class:`dict`.
    	"""
        if not excluded:
            excluded = []
        color_groups = collections.defaultdict(list)
        for atom in self.atoms:
            if atom.atom_number not in excluded:
                color_groups[atom.color].append(atom.atom_number)
        return color_groups

    def detect_abnormal_atom(self):
        """
        To find the atoms with invalid bond counts.

        :return: a list of atom numbers with invalid bond counts.
        :rtype: :py:class:`list`.
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

    def curate_invalid_n(self):
        """
        To curate charge of the invalid N atoms.

        :return: None.
        :rtype: :py:obj:`None`.
    	"""
        abnormal_atoms = self.detect_abnormal_atom()
        for atom_index in abnormal_atoms["N"]:
            atom = self.atoms[atom_index]
            atom.charge += 1
            for neighbor in atom.neighbors:
                if self.atoms[neighbor].default_symbol == "O" and self.bond_lookup[(atom.atom_number, neighbor)].bond_type == "2":
                    self.atoms[neighbor].charge -= 1
                    self.bond_lookup[(atom.atom_number, neighbor)].update_stereochemistry("1")

    # def metals(self):
    #     """
    #     To get all the metals in the compound.
    #
    #     :return: the dictionary of metal and the corresponding list of index.
    # 	"""
    #     metals = collections.defaultdict(list)
    #     for index, atom in enumerate(self.atoms):
    #         if atom.default_symbol in metal_symbols:
    #             metals[atom.default_symbol].append(index)
    #     return metals

    # def update_aromatic_bond_type(self, cycles):
    #     """
    #     Update the aromatic bond types.
    #     Two cases: 1) change the bond in the aromatic ring to aromatic bond (bond type = 4)
    #            2) change the double bond connecting to the aromatic ring to single bond.
    #     :param list cycles: the list of cycles of aromatic atom index.
    #
    #     :return:
    #     """
    #     atom_in_cycle = [atom for cycle in cycles for atom in cycle]
    #     for cycle in cycles:
    #         aromatic_bonds = self.extract_aromatic_bonds(cycle)
    #         for bond in aromatic_bonds:
    #             bond.update_bond_type("4")
    #     bond_out_of_cycle = self.extract_double_bond_connecting_cycle(atom_in_cycle)
    #     for bond in bond_out_of_cycle:
    #         bond.update_bond_type("1")

    def update_aromatic_bond(self, aromatic_bonds, aromatic_atoms):
        """
        To update the bond type of aromatic bonds in the compound.
        Two steps are involved:
            1) Update the bond in the aromatic ring to type "4".
            2) Update the double bond connecting to the aromatic ring to "1".

        :param aromatic_bonds: a list of aromatic bonds in the compound (represented by the first and second atom number).
        :type aromatic_bonds: :py:class:`list`.
        :param aromatic_atoms: a list of aromatic atoms in the aromatic ring.
        :type aromatic_atoms: :py:class:`list`.
        :return: None.
        :rtype: :py:obj:`None`.
        """

        for atom_i, atom_j in aromatic_bonds:
            self.bond_lookup[(atom_i, atom_j)].update_bond_type("4")
        bond_out_of_cycle = self.extract_double_bond_connecting_cycle(aromatic_atoms)
        for bond in bond_out_of_cycle:
            bond.update_bond_type("1")

    def extract_double_bond_connecting_cycle(self, atom_in_cycle):
        """
        To extract the double bonds connecting to the aromatic cycles.

        :param atom_in_cycle: the list of atoms in the aromatic cycles.
        :type atom_in_cycle: :py:class:`list`.
        :return: the list of outside double bond connecting to the aromatic cycles.
        :rtype: :py:class:`list`.
        """
        double_bond_connecting_cycle = []
        for atom_index in atom_in_cycle:
            if self.atoms[atom_index].default_symbol == "C":
                for neighbor_index in self.atoms[atom_index].neighbors:
                    if neighbor_index not in atom_in_cycle and self.bond_lookup[(atom_index, neighbor_index)].bond_type == "2":
                        double_bond_connecting_cycle.append(self.bond_lookup[(atom_index, neighbor_index)])
        return double_bond_connecting_cycle

    # def extract_aromatic_bonds(self, cycle):
    #     """
    #     Extract the aromatic bonds based on the atom indexes in the cycle.
    #     :param list cycle: the
    #
    #     :return: the list of aromatic bond.
    # 	"""
    #     aromatic_bonds = []
    #     all_pairs = list(itertools.combinations(cycle, 2))
    #     visited = set()
    #     for pair in all_pairs:
    #         if (pair[0], pair[1]) not in visited and (pair[0], pair[1]) in self.bond_lookup:
    #             aromatic_bonds.append(self.bond_lookup[(pair[0], pair[1])])
    #             visited.add((pair[0], pair[1]))
    #             visited.add((pair[1], pair[1]))
    #     return aromatic_bonds

    def separate_connected_components(self, bonds):
        """
        This is used in constructing the aromatic substructures detected by Indigo method.
        A compound can have several disjoint aromatic substructures. Here, we need to find the disjoint parts.
        The basic idea is union-find. We union atoms that are connected by a bond.

        :param bonds: the list of bonds representing by the atom numbers in the bond.
        :type bonds: :py:class:`list`.
        :return: a list of components represented by a list atom numbers in each component.
        :rtype: :py:class:`list`.
        """

        atoms = set()
        for atom_i, atom_j in bonds:
            atoms.add(atom_i)
            atoms.add(atom_j)

        parent_index = {i: i for i in atoms}

        def find_parent(i):
            if i != parent_index[i]:
                parent_index[i] = find_parent(parent_index[i])
            return parent_index[i]

        def union(p_1, p_2):
            parent_index[p_1] = p_2

        for atom_i, atom_j in bonds:
            p_1 = find_parent(atom_i)
            p_2 = find_parent(atom_j)
            if p_1 != p_2:
                union(p_1, p_2)

        groups = collections.defaultdict(list)
        for i in atoms:
            p_i = find_parent(i)
            groups[p_i].append(i)
        return [groups[key] for key in groups]

    def connected_components(self):
        """
        Detect the connected components in the compound structure. (using the breadth first search)
        Cases when not all the atoms are connected together.

        :return: the dictionary of the connected components.
        :rtype: :py:class:`dict`.
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

    def calculate_distance_to_r_groups(self):
        """
        To calculate the distance of each atom to its nearest R group (using the dijkstra's algorithm).

        :return: None:
        :rtype: :py:obj:`None`.
        """
        distance_matrix = [len(self.heavy_atoms)] * len(self.heavy_atoms) 

        if self.r_groups:
            for r_index in self.r_groups:
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
                self.heavy_atoms[i].distance_to_r = dist

    def find_cycles(self, short_circuit=False, cutoff=40):
        """
        Find the cycles in the compound.

        :return: the list of cycles in the compound.
        :rtype: :py:class:`list`.
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
        return [list(x) for x in set(tuple(x) for x in [sorted(i[:-1]) for l in all_cycles for i in l])]

    def structure_matrix(self, resonance=False, backbone=False):
        """
         To construct graph structural matrix of this compound.
         matrix[i][j] = 0 suggests the two atoms are not connected directly.
         Other integer represented the bond type connecting the two atoms.

        :param resonance: bool whether to ignore the difference between single and double bonds.
        :type resonance: :py:obj:`bool`.
        :param backbone: bool whether to ignore bond types. This is for parsing atoms mappings from KEGG RCLASS.
        :type backbone: :py:obj:`bool`.
        :return: the constructed structure matrix for this compound.
        :rtype: :py:class:`ndarray`.
        """
        matrix = numpy.zeros((len(self.heavy_atoms), len(self.heavy_atoms)), dtype=numpy.uint8)
        for bond in self.bonds:
            atom_1, atom_2 = bond.first_atom_number, bond.second_atom_number
            if atom_1 in self.index_of_heavy_atoms and atom_2 in self.index_of_heavy_atoms:
                bond_type = int(bond.bond_type)
                if resonance and bond_type == 2:
                    bond_type = 1
                if backbone:
                    bond_type = 1
                matrix[self.index_of_heavy_atoms[atom_1]][self.index_of_heavy_atoms[atom_2]] = bond_type
                matrix[self.index_of_heavy_atoms[atom_2]][self.index_of_heavy_atoms[atom_1]] = bond_type
        return matrix

    @property
    def distance_matrix(self):
        """
        To construct the distance matrix of the compound. (using the Floyd Warshall Algorithm)
        distance[i][j] suggests the distance between atom i and j.

        :return: the distance matrix of the compound.
        :rtype: :py:class:`ndarray`.
    	"""
        if self.heavy_atoms:
            distance_matrix = numpy.ones((len(self.heavy_atoms), len(self.heavy_atoms)),dtype=numpy.uint16 ) * len(self.heavy_atoms)
            for bond in self.bonds:
                atom_1, atom_2 = bond.first_atom_number, bond.second_atom_number
                if atom_1 in self.index_of_heavy_atoms and atom_2 in self.index_of_heavy_atoms:
                    distance_matrix[self.index_of_heavy_atoms[atom_1]][self.index_of_heavy_atoms[atom_2]] = 1
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
        To update the color tuple of the atoms in the compound. This color tuple includes information of its neighboring
        atoms and bonds.

        :param resonance: bool whether to ignore the difference between single and double bonds.
        :type resonance: :py:obj:`bool`.
        :return: None.
        :rtype: :py:obj:`None`.
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
                atom.color_tuple = (elements["C"], elements["N"], elements["S"], elements["O"], bond_types["1"] +
                                    bond_types["2"])
            else:
                atom.color_tuple = (elements["C"], elements["N"], elements["S"], elements["O"], bond_types["1"],
                                    bond_types["2"])
    
    def find_mappings(self, the_other, resonance=True, r_distance=False, backbone=False):
        """
        Find the one to one atom mappings between two compounds using BASS algorithm.

        :param the_other: the mappings compound entity.
        :type the_other: :class:`~MDH.compound.Compound`.
        :param resonance: whether to ignore the difference between single and double bonds.
        :type resonance: :py:obj:`bool`.
        :param r_distance: whether to take account of position of R groups.
        :type r_distance: :py:obj:`bool`.
        :param backbone: whether to ignore the bond types.
        :type backbone: :py:obj:`bool`.
        :return: the list of atom mappings in the heavy atom order.
        :rtype: :py:class:`list`.
    	"""
        self.update_color_tuple(resonance=resonance)
        the_other.update_color_tuple(resonance=resonance)
        mappings = []
        mapping_matrix = BASS.make_mapping_matrix(the_other, self, True, True, r_distance)
    
        if mapping_matrix is not None:
            mappings = BASS.find_mappings(the_other.structure_matrix(resonance=resonance, backbone=backbone), the_other.distance_matrix,
                                          self.structure_matrix(resonance=resonance, backbone=backbone), self.distance_matrix, mapping_matrix)
        # for the mappings, the from_idx, to_idx in enumerate(mapping), from_idx is in the_other_compound, to_idx is in the self.
        one_to_one_mappings = []
        for sub in mappings:
            cur_mappings = {}
            for from_index, to_index in enumerate(sub):
                cur_mappings[self.heavy_atoms[to_index].atom_number] = the_other.heavy_atoms[from_index].atom_number
            one_to_one_mappings.append(cur_mappings)
        return one_to_one_mappings

    def map_resonance(self, the_other, r_distance=False):
        """
        Check if the resonant mappings are valid between the two compound structures. If the mapped atoms don't share
        the same local coloring identifier, we check if the difference is caused by the position of double bonds.
        Find the three atoms involved in the resonant structure and check if one of the atom is not C.
                N (a)            N (a)
                / \\             // \
            (b) C   N (c)    (b) C   N (c)

        :param the_other: the mappings compound entity.
        :type the_other: :class:`~MDH.compound.Compound`.
        :param r_distance: to take account of positon of R groups.
        :type r_distance: :py:obj:`bool`.
	    :return: the list of valid atom mappings between the two compound structures.
	    :rtype: :py:class:`list`.
    	"""
        one_to_one_mappings = self.find_mappings(the_other, resonance=True, r_distance=r_distance)
        valid_mappings = []
        for cur_mappings in one_to_one_mappings:
            flag = True
            reversed_mappings = { cur_mappings[key] : key for key in cur_mappings }
            for from_index in sorted(cur_mappings.keys()):
                to_index = cur_mappings[from_index]

                # Detect if the two atoms have the same local binding environment by comparing the first layer color identifier.
                if 1 in self.atoms[from_index].color_layers and 1 in the_other.atoms[to_index].color_layers and \
                        self.atoms[from_index].color_layers[1] == the_other.atoms[to_index].color_layers[1]:
                    continue
                
                three = {from_index}
                from_dbond_atom_index = self.find_double_bond_linked_atom(from_index)
                to_dbond_atom_index = the_other.find_double_bond_linked_atom(to_index)
                
                #cannot find directly linked double bond.
                if from_dbond_atom_index == -1 and to_dbond_atom_index == -1:
                    flag = False
                #Need to conduct second search
                # atom b in the picture
                elif from_dbond_atom_index == -1:
                    reversed_from_dbond_atom_index = reversed_mappings[to_dbond_atom_index] # here we can find the atom a
                    three.add(reversed_from_dbond_atom_index)
                    from_next_dbond_atom_index = self.find_double_bond_linked_atom(reversed_from_dbond_atom_index) 
                    # based on atom a we can find atom c
                    if from_next_dbond_atom_index != -1:
                        three.add(from_next_dbond_atom_index)

                # atom c in the picture
                elif to_dbond_atom_index == -1 and from_dbond_atom_index in cur_mappings:
                    three.add(from_dbond_atom_index) # here we find atom a
                    reversed_to_dbond_atom_index = cur_mappings[from_dbond_atom_index]
                    # based on atom a find atom b in the other compound.
                    to_next_dbond_atom_index = the_other.find_double_bond_linked_atom(reversed_to_dbond_atom_index)
                    if to_next_dbond_atom_index != -1:
                        three.add(reversed_mappings[to_next_dbond_atom_index])

                # atom a in the picture
                # find the other two atoms directly.
                elif from_dbond_atom_index != -1 and to_dbond_atom_index != -1:
                    three.add(from_dbond_atom_index)
                    three.add(reversed_mappings[to_dbond_atom_index])

                if len(three) < 3:
                    flag = False
                else:
                    k = [self.atoms[i].default_symbol for i in three].count("C")
                    if k == 3:
                        flag = False

                if not flag:
                    break
            if flag:
                valid_mappings.append(cur_mappings)
        return valid_mappings

    def find_double_bond_linked_atom(self, i):
        """
        Find the atom that is doubly linked to the target atom[i].

        :param i: the ith atom in the compound.
        :type i: :py:class:`int`.
        :return: the index of the doubly linked atom.
        :rtype: :py:class:`int`.
        """
        for neighbor_index in self.atoms[i].neighbors:
            if neighbor_index in self.index_of_heavy_atoms and self.bond_lookup[(i, neighbor_index)].bond_type == "2":
                return neighbor_index
        return -1

    def define_bond_stereochemistry(self):
        """
        Define the stereochemistry of double bonds in the compound.

        :return: None.
        :rtype: :py:obj:`None`.
    	"""
        for bond in self.bonds:
            if bond.bond_type == "2":
                bond.update_stereochemistry(self.calculate_bond_stereochemistry(bond))
    
    def calculate_bond_stereochemistry(self, bond):
        """
        Calculate the stereochemisty of double bond based on the geometric properties. The line of double bond divides
        the plane into two parts. For the atom forming the doble bond, it normally has two branches. If the two branches
        are not the same, we can them heavy side and light side (heavy side containing atoms with heavier atomic weights).
        We determine the bond stereochemistry by checking if the two heavy sides lie on the same part of the divided plane.

            H   L   H   H
            \___/   \___/
             ___     ___
            /   \   /   \
            L   H   L   L
            trans   cis
        :param bond: the bond entity.
        :type bond: :class:`~MDH.compound.Bond`.
        :return: the calculated bond stereochemistry.
        :rtype: :py:class:`int`.
    	"""
        vertical = False

        first_atom = self.atoms[bond.first_atom_number]
        second_atom = self.atoms[bond.second_atom_number]
        slope, b = 0, 0
        if first_atom.x != second_atom.x:
            slope = (first_atom.y-second_atom.y)/(first_atom.x-second_atom.x)
            b = first_atom.y - slope * first_atom.x 
        else:
            vertical = True

        first_atom_neighbor_index = [neighbor for neighbor in first_atom.neighbors if neighbor != second_atom.atom_number]
        second_atom_neighbor_index = [neighbor for neighbor in second_atom.neighbors if neighbor != first_atom.atom_number]
        
        if not first_atom_neighbor_index or not second_atom_neighbor_index:
            return 0
        
        first_heavy_branch, first_light_branch = self.compare_branch_weights(first_atom_neighbor_index, first_atom)
        second_heavy_branch, second_light_branch = self.compare_branch_weights(second_atom_neighbor_index, second_atom)

        if (first_light_branch and first_light_branch == first_heavy_branch) or (second_light_branch and second_light_branch == second_heavy_branch):
            return 0
       
        if vertical:
            if (self.atoms[first_heavy_branch].x - first_atom.x) * (self.atoms[second_heavy_branch].x - second_atom.x) > 0:
                return 1
            else:
                return -1
        else:
            if (self.calculate_y_coordinate(slope, b, self.atoms[first_heavy_branch]) - self.atoms[first_heavy_branch].y) *\
                    (self.calculate_y_coordinate(slope, b, self.atoms[second_heavy_branch]) - self.atoms[second_heavy_branch].y) > 0:
                return 1
            else:
                return -1

    @staticmethod
    def calculate_y_coordinate(slope, b, atom):
        """
        Calculate the y coordinate of the atom based on the linear function: y = slope * x + b

        :param slope: the slope of the targeted line.
        :type slope: :py:class:`float`.
        :param b: the intercept of the targeted line.
        :type b: :py:class:`float`.
        :param atom: the atom entity.
        :type atom: :class:`~MDH.compound.Atom`
        :return: the calculated y coordinate.
        :rtype: :py:class:`float`.
    	"""
        return atom.x * slope + b

    def collect_atomic_weights_of_neighbors(self, neighbors):
        """
        To collect the atomic weights of the current layer's neighbors.

        :param neighbors: the list of neighboring atom numbers.
        :type neighbors: :py:class:`list`.
        :return: the list of atomic weight for this layer's neighbors.
        :rtype: :py:class:`list`.
    	"""
        neighbor_atomic_weights = [atomic_weights[self.atoms[index].default_symbol] for index in neighbors]
        neighbor_atomic_weights.sort(reverse=True)
        return neighbor_atomic_weights

    def compare_branch_weights(self, neighbors, atom_forming_double_bond):
        """
        To determine the heavy and light branches connecting the atom forming the double bond. This is based on
        comparison of the atomic weights of the two branches. (Breadth first algorithm).

        :param neighbors: the list of atom numbers of the atoms connecting the atom forming the double bond.
        :type neighbors: :py:class:`list`.
        :param atom_forming_double_bond: the bond entity.
	    :type atom_forming_double_bond: :class:`~MDH.compound.Bond`
        :return: heavy and light branches. [heavy_side, light_side]
        :rtype: :py:class:`list`.
        """
        if len(neighbors) < 2:
            return neighbors[0], None

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

    def get_next_layer_neighbors(self, cur_layer_neighbors, visited, excluded=None):
        """
        To get the next layer's neighbors.

        :param cur_layer_neighbors: the list of atom numbers of the current layer.
        :type cur_layer_neighbors: :py:class:`list`.
        :param visited: the atom numbers that have already visited.
        :type visited: :py:class:`set`.
        :param excluded: the list of atom numbers that should not be included in the next layer.
        :type excluded: :py:class:`list`.
        :return: the atom numbers of next layer's neighbors.
        :rtype: :py:class:`list`.
        """
        if not excluded:
            excluded = []
        next_layer_neighbors = []
        for index in cur_layer_neighbors:
            for next_neighbor in self.atoms[index].neighbors:
                if next_neighbor not in visited and next_neighbor not in excluded:
                    next_layer_neighbors.append(next_neighbor)
                    visited.add(next_neighbor)
        return next_layer_neighbors

    def color_compound(self, r_groups=True, bond_stereo=False, atom_stereo=False, resonance=False, isotope_resolved=False,
                       charge=False, backbone=False):
        """
        To color the compound.

        :param r_groups:  If true, add R groups in the coloring.
        :type r_groups: :py:obj:`bool`.
        :param bond_stereo:  If true, add stereo information to bonds when constructing colors.
        :type bond_stereo: :py:obj:`bool`.
        :param atom_stereo: If true, add atom stereo information when constructing colors.
        :type atom_stereo: :py:obj:`bool`.
        :param resonance: If true, ignore difference between double bonds and single bonds.
        :type resonance: :py:obj:`bool`.
        :param isotope_resolved: If true, add isotope information when constructing colors.
        :type isotope_resolved: :py:obj:`bool`.
        :param charge: If true, add charge information when constructing colors.
        :type charge: :py:obj:`bool`.
        :param backbone: If true, ignore bond types in the coloring.
        :type backbone: :py:obj:`bool`.
        :return: None.
        :rtype: :py:obj:`None`.
    	"""
        excluded_index = self.metal_index + self.h_index
        excluded_index += self.r_groups if r_groups else []
        atoms_to_color = [i for i in range(len(self.atoms)) if i not in excluded_index]
        self.reset_color()
        self.generate_atom_zero_layer_color(isotope_resolved=isotope_resolved, charge=charge, atom_stereo=atom_stereo)
        self.first_round_color(atoms_to_color, excluded_index=excluded_index, bond_stereo=bond_stereo, resonance=resonance,
                               backbone=backbone)
        self.curate_invalid_symmetric_atoms(atoms_to_color, excluded_index=excluded_index, bond_stereo=bond_stereo,
                                            resonance=resonance, backbone=backbone)
        self.color_metal(bond_stereo=bond_stereo, resonance=resonance)
    
    def reset_color(self):
        """
        To set the color of atoms and compound to empty.

        :return: None:
        :rtype: :py:obj:`None`.
    	"""
        for atom in self.atoms:
            atom.reset_color()

    def generate_atom_zero_layer_color(self, isotope_resolved=False, charge=False, atom_stereo=False):
        """
        To generate the zero layer color identifier for each atom. We don't consider H and metals here.

        :param isotope_resolved: If true, add isotope information when constructing colors.
        :type isotope_resolved: :py:obj:`bool`.
        :param charge: If true, add charge information when constructing colors.
        :type charge: :py:obj:`bool`.
        :param atom_stereo: If true, add atom stereochemistry information when constructing colors.
        :type atom_stereo: :py:obj:`bool`.
        :return: None.
        :rtype: :py:obj:`None`.
    	"""
        for index, atom in enumerate(self.atoms):
            atom.color_atom(isotope_resolved=isotope_resolved, charge=charge, atom_stereo=atom_stereo)

    def generate_atom_color_with_neighbors(self, atom_index, excluded=None, zero_core_color=True, zero_neighbor_color=True,
                                           resonance=False, bond_stereo=False, backbone=False):
        """
        To generate the atom color with its neighbors. We add this color name when we try to incorporate neighbors'
        information in naming.
        Here, we don't need to care about the atom stereo. It has been taken care of in generating color_0.
        Basic color formula: atom.color + [neighbor.color + bond.bond_type]

        :param list atom_index: the list of atom index to color.
        :param list excluded: the list of atom index will be excluded from coloring.
        :param zero_core_color: If ture, we use the atom.color_0 else atom.color for the core atom.
        :type zero_core_color: :py:obj:`bool`.
        :param zero_neighbor_color: If ture, we use the atom.color_0 else atom.color for the neighbor atoms.
        :type zero_neighbor_color: :py:obj:`bool`.
        :param resonance: If true, detect resonant compound pairs without distinguish between double bonds and single bonds.
        :type resonance: :py:obj:`bool`.
        :param bond_stereo:  If true, add stereo information to bonds when constructing colors.
        :type bond_stereo: :py:obj:`bool`.
        :param backbone: If true, ignore bond types in the coloring.
        :type backbone: :py:obj:`bool`.
        :return: the dictionary of atom index to its color name containing neighbors.
        :rtype: :py:class:`dict`.
    	"""
        if not excluded:
            excluded = []
        atom_color_with_neighbors = collections.defaultdict(str)
        for index in atom_index:
            atom = self.atoms[index]
            atom_color = atom.color_0 if zero_core_color else atom.color
            color_elements = collections.Counter()
            for neighbor_index in atom.neighbors:
                if neighbor_index not in excluded:
                    neighbor_color = self.atoms[neighbor_index].color_0 if zero_neighbor_color else self.atoms[neighbor_index].color
                    connecting_bond = self.bond_lookup[(atom.atom_number, neighbor_index)]
                    bond_type = connecting_bond.bond_type
                    if resonance and bond_type == "2":
                        bond_type = "1"
                    if backbone:
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

    def first_round_color(self, atoms_to_color, excluded_index=None, bond_stereo=False, resonance=False, backbone=False, depth=5000):
        """
        To do the first round of coloring this compound. We add neighbors' information layer by layer to the atom color
        identifier until it has a unique identifier or all the atoms in the compound have been used for naming.
        (based on the breadth first search algorithm)

        :param atoms_to_color: the list of atom numbers of atoms to be colored.
        :type atoms_to_color: :py:class:`list`.
        :param excluded_index: the list of atom numbers of atoms to be excluded from coloring.
        :type excluded_index: :py:class:`list`.
        :param bond_stereo: If true, add stereo information to bonds when constructing colors.
        :type bond_stereo: :py:obj:`bool`.
        :param resonance: If true, ignore the difference between double bonds and single bonds.
        :type resonance: :py:obj:`bool`.
        :param backbone: If true, ignore bond types in the coloring.
        :type backbone: :py:obj:`bool`.
        :param depth: the max depth of coloring.
        :type depth: :py:class:`int`.
        :return: None.
        :rtype: :py:obj:`None`.
        """
        atom_color_with_neighbors = self.generate_atom_color_with_neighbors(atoms_to_color, excluded=excluded_index,
                                                                            zero_core_color=True, zero_neighbor_color=True,
                                                                            resonance=resonance, bond_stereo=bond_stereo,
                                                                            backbone=backbone)
        if not excluded_index:
            excluded_index = []

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
                atom = self.atoms[atom_index]
                color_elements = collections.Counter()
                for neighbor_index in atom_neighbors[atom_index]:
                    if neighbor_index not in excluded_index:
                        color_elements[atom_color_with_neighbors[neighbor_index]] += 1
                added = ""
                for name in sorted(color_elements):
                    added += "({0}_{1})".format(color_elements[name], name) if color_elements[name] > 1 else "({0})".format(name)
                atom.color += added
                atom.color_layers[i+1] = atom.color
                atom_neighbors[atom.atom_number] = self.get_next_layer_neighbors(atom_neighbors[atom_index], visited[atom_index],
                                                                                 excluded=excluded_index)
                current_layer_color_groups[atom.color].append(atom_index)

            if i > 3:
                # avoid early stop
                atom_to_color_update = []
                for name in current_layer_color_groups.keys():
                    if len(current_layer_color_groups[name]) > 1:
                        atom_to_color_update.extend(current_layer_color_groups[name])
                atoms_to_color = atom_to_color_update
            i += 1
    
    def invalid_symmetric_atoms(self, atoms_to_color, excluded_index=None, bond_stereo=False, resonance=False, backbone=False):
        """
        To check if the atoms with the same color identifier are symmetric.

        :param atoms_to_color: the list of atom numbers of atoms to be colored.
        :type atoms_to_color: :py:class:`list`.
        :param excluded_index: the list of atom numbers of atoms to be excluded from coloring.
        :type excluded_index: :py:class:`list`.
        :param bond_stereo: If true, add stereo information to bonds when constructing colors.
        :type bond_stereo: :py:obj:`bool`.
        :param resonance: If true, ignore the difference between double bonds and single bonds.
        :type resonance: :py:obj:`bool`.
        :param backbone: If true, ignore bond types in the coloring.
        :type backbone: :py:obj:`bool`.
        :return: the list of atoms to be recolored.
        :rtype: :py:class:`list`.
        """
        atom_color_with_neighbors = self.generate_atom_color_with_neighbors(atoms_to_color, excluded=excluded_index,
                                                                            zero_core_color=False, zero_neighbor_color=False,
                                                                            resonance=resonance, bond_stereo=bond_stereo,
                                                                            backbone=backbone)
        not_valid = []
        color_groups = self.color_groups(excluded=excluded_index)
        if not excluded_index:
            excluded_index = []
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
                    target_color_list = [atom_color_with_neighbors[atom_index] for atom_index in current_layer[target_index]
                                         if atom_index not in excluded_index]
                    target_color_list.sort()
                    target_color = ''.join(target_color_list)

                    #check if neighbors of this layer are the same among all the atoms with the same color identifier.
                    for compared_index in current_layer.keys():
                        if compared_index != target_index:
                            compared_color_list = [atom_color_with_neighbors[atom_index] for atom_index in
                                                   current_layer[compared_index] if atom_index not in excluded_index]
                            compared_color_list.sort()
                            compared_atom_color = ''.join(compared_color_list)

                            if target_color != compared_atom_color:
                                not_valid.append(atom_index_with_same_color)
                                flag = False
                                break
                    # If the share the same neighbors, check the next layer.
                    if flag:
                        for atom_index in current_layer.keys():
                            next_layer_neighbors = self.get_next_layer_neighbors(current_layer[atom_index], visited[atom_index],
                                                                                 excluded=excluded_index)
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

    def curate_invalid_symmetric_atoms(self, atoms_to_color, excluded_index=None, bond_stereo=False, resonance=False, backbone=False):
        """
        To curate the atom color identifier of invalid symmetric atom.
        We recolor those invalid atoms with the full color identifiers of its neighbors layer by layer until where the
        difference can be captured.

        :param atoms_to_color: the list of atom numbers of atoms to be colored.
        :type atoms_to_color: :py:class:`list`.
        :param excluded_index: the list of atom numbers of atoms to be excluded from coloring.
        :type excluded_index: :py:class:`list`.
        :param bond_stereo: If true, add stereo information to bonds when constructing colors.
        :type bond_stereo: :py:obj:`bool`.
        :param resonance: If true, ignore the difference between double bonds and single bonds.
        :type resonance: :py:obj:`bool`.
        :param backbone: If true, ignore bond types in the coloring.
        :type backbone: :py:obj:`bool`.
        :return: None.
        :rtype: :py:obj:`None`.
        """
        not_valid = self.invalid_symmetric_atoms(atoms_to_color, excluded_index, bond_stereo=bond_stereo, resonance=resonance)
        while not_valid:
            atom_color_with_neighbors = self.generate_atom_color_with_neighbors(atoms_to_color, excluded=excluded_index,
                                                                                zero_core_color=False, zero_neighbor_color=False,
                                                                                resonance=resonance, bond_stereo=bond_stereo,
                                                                                backbone=backbone)
            for invalid_symmetric_atom_index in not_valid:
                for atom_index in invalid_symmetric_atom_index:
                    visited = {atom_index}
                    current_layer_neighbors = [atom_index]
                    atom_color = ""
                    while current_layer_neighbors:
                        color_elements = collections.Counter()
                        for neighbor_index in current_layer_neighbors:
                            if neighbor_index not in excluded_index:
                                color_elements[atom_color_with_neighbors[neighbor_index]] += 1
                        added = ""
                        for name in sorted(color_elements):
                            added += "({0}_{1})".format(color_elements[name], name) if color_elements[name] > 1 else \
                                "({0})".format(name)
                        atom_color += added
                        current_layer_neighbors = self.get_next_layer_neighbors(current_layer_neighbors, visited,
                                                                                excluded=excluded_index)
                    self.atoms[atom_index].color = atom_color
            not_valid = self.invalid_symmetric_atoms(atoms_to_color, excluded_index, bond_stereo=bond_stereo,
                                                     resonance=resonance)

    def color_metal(self, bond_stereo=False, resonance=True, backbone=False):
        """
        To color the metals in the compound. Here we just incorporate information of directly connected atoms.

        :param bond_stereo: If true, add stereo information to bonds when constructing colors.
        :type bond_stereo: :py:obj:`bool`.
        :param resonance: If true, ignore difference between double bonds and single bonds.
        :type resonance: :py:obj:`bool`.
        :param backbone: If true, ignore bond types.
        :type backbone: :py:obj:`bool`.
        :return: None.
        :rtype: :py:obj:`None`.
    	"""
        atom_color_with_neighbors = self.generate_atom_color_with_neighbors(self.metal_index, excluded=self.h_index,
                                                                            zero_core_color=True, zero_neighbor_color=False,
                                                                            resonance=resonance, bond_stereo=bond_stereo,
                                                                            backbone=backbone)
        for atom_index in self.metal_index:
            self.atoms[atom_index].color = atom_color_with_neighbors[atom_index]
      
    def color_h(self, bond_stereo=False, resonance=True, backbone=False):
        """To color the metals in the compound. Here we just incorporate information of directly connected atoms.

        :param bond_stereo:  If true, add stereo information to bonds when constructing colors.
        :type bond_stereo: :py:obj:`bool`.
        :param resonance: If true, detect resonant compound pairs without distinguish between double bonds and single bonds.
        :type resonance: :py:obj:`bool`.
        :param backbone: If true, ignore bond types.
        :type backbone: :py:obj:`bool`.
        :return: None.
        :rtype: :py:obj:`None`.
    	"""
        atom_color_with_neighbors = self.generate_atom_color_with_neighbors(self.h_index, zero_core_color=True,
                                                                            zero_neighbor_color=False, resonance=resonance,
                                                                            bond_stereo=bond_stereo, backbone=backbone)
        for atom_index in self.h_index:
            self.atoms[atom_index].color = atom_color_with_neighbors[atom_index]

    def metal_color_identifier(self, details=True):
        """
        To generate the metal color string representation.

        :param details: if true, use full metal color when constructing identifier.
        :type details: :py:obj:`bool`.
        :return: the metal color string representation.
        :rtype: :py:class:`str`.
    	"""
        if details:  
            color_counter = collections.Counter([self.atoms[index].color for index in self.metal_index])
        else:
            color_counter = collections.Counter([self.atoms[index].color_0 for index in self.metal_index])
        return "".join(["({0})({1})".format(color_counter[key], key) for key in sorted(color_counter)])
    
    def h_color_identifier(self, details=True):
        """
        To generate the H color string representation.

        :param details: if true, use full H color when constructing identifier.
        :type details: :py:obj:`bool`.
        :return: the H color string representation.
        :rtype: :py:class:`str`.
    	"""
        if details:
            color_counter = collections.Counter([self.atoms[index].color for index in self.h_index])
        else:
            color_counter = collections.Counter([self.atoms[index].color_0 for index in self.h_index])
        return "".join(["({0})({1})".format(color_counter[key], key) for key in sorted(color_counter)])

    def backbone_color_identifier(self, r_groups=False):
        """
        To generate the backbone color identifier for this compound. Exclude Hs and metals.

        :param r_groups: bool whether to include the R group.
        :type r_groups: :py:obj:`bool`.
        :return: the color string identifier for this compound.
        :rtype: :py:class:`str`.
        """
        excluded_index = self.metal_index + self.h_index
        if r_groups:
            excluded_index += self.r_groups
        color_groups = self.color_groups(excluded=excluded_index)
        return "".join(["({0})({1})".format(len(color_groups[key]), key) for key in sorted(color_groups)])

    # To detect if this compound can be paired to the other compound.
    # Three cases: 1) Both compounds are specific compounds (no R groups); Just compare the coloring identifier.
    #              2) One compound contains R group(s);
    #              3) Both compounds contain R group(s).
    # If the two compounds can be paired, we need to determine their relationship by checking the chemical details (eg:
    # bond stereochemistry and atom stereochemistry.


    def get_chemical_details(self, excluded=None):
        """
        To get the chemical details of the compound, which include the atom stereo chemistry and bond stereo chemistry.
        This is to compare compound with the same structures (or the same color identifiers).

        :param excluded: the atoms to be ignored.
        :type excluded: :py:class:`list`.
        :return: the list of chemical details in the compound.
        :rtype: :py:class:`list`.
        """
        if not excluded:
            excluded = []
        chemical_details = []
        for atom in self.atoms:
            if atom.atom_number not in excluded and atom.atom_stereo_parity != "0":
                chemical_details.append("{0}-{1}".format(atom.color, atom.atom_stereo_parity))
        for bond in self.bonds:
            if bond.first_atom_number not in excluded and bond.second_atom_number not in excluded and \
                    bond.bond_stereo != "0":
                names = [self.atoms[bond.first_atom_number].color, self.atoms[bond.second_atom_number].color]
                names.sort()
                chemical_details.append("{0}-{1}-{2}".format(names[0], names[1], bond.bond_stereo))
        chemical_details.sort()
        return chemical_details

    @staticmethod
    def compare_chemical_details(one_chemical_details, the_other_chemical_details):
        """
        To compare the chemical details of the two compounds.
        Then return the relationship between the two compounds.
        The relationship can be equivalent, generic-specific and loose, represented by 0, (-1, 1), 2

        :param one_chemical_details: the chemical details of one compound.
        :type one_chemical_details: :py:class:`list`.
        :param the_other_chemical_details: the chemical details of the other compound.
        :type the_other_chemical_details: :py:class:`list`.
        :return: the relationship between the two structures and the count of not mapped chemical details.
        :rtype: :py:class:`int`.
        """
        # order is small to big
        one_more = []
        the_other_more = []
        while one_chemical_details and the_other_chemical_details:
            if one_chemical_details[-1] == the_other_chemical_details[-1]:
                one_chemical_details.pop()
                the_other_chemical_details.pop()
            elif one_chemical_details[-1] > the_other_chemical_details[-1]:
                one_more.append(one_chemical_details.pop())
            else:
                the_other_more.append(the_other_chemical_details.pop())
        one_more.extend(one_chemical_details)
        the_other_more.extend(the_other_chemical_details)
        if not one_more and not the_other_more:
           return 0, 0
        elif one_more and the_other_more:
            return 2, len(one_more) + len(the_other_more)
        elif one_more:
            return 1, len(one_more)# one_compound is more specific than the_other_compound
        elif the_other_more:
            return -1, len(the_other_more) # the_other_compound is more specific than one_compound

    def same_structure_relationship(self, the_other_compound):
        """
        To determine the relationship of two compounds with the same structure.

        :param the_other_compound: the other :class:`~MDH.compound.Compound` entity.
        :type the_other_compound: :class:`~MDH.compound.Compound`.
        :return: the relationship and the atom mappings between the two compounds.
        :rtype: :py:class:`int` and :py:class:`dict`.
        """
        relationship, unmapped_count = self.compare_chemical_details(self.get_chemical_details(),
                                                                the_other_compound.get_chemical_details())
        # return relationship and atom mappings.
        return relationship, self.generate_atom_mapping_by_atom_color(the_other_compound)

    def generate_atom_mapping_by_atom_color(self, the_other_compound):
        """
        To generate the atom mappings between the two compounds.
        Assume the two compounds have the same structure, so we can achieve atom mappings through atom colors.

        :param the_other_compound: the other :class:`~MDH.compound.Compound` entity.
        :type the_other_compound: :class:`~MDH.compound.Compound`.
        :return: the atom mappings between the two compounds.
        :rtype: :py:class:`dict`.
        """

        one_to_one_mappings = collections.defaultdict(list)
        for idx, atom in enumerate(self.atoms):
            for the_other_idx, the_other_atom in enumerate(the_other_compound.atoms):
                if atom.color == the_other_atom.color:
                    one_to_one_mappings[idx].append(the_other_idx)
        return one_to_one_mappings

    def optimal_resonant_mapping(self, the_other_compound, mappings):
        """

        :param the_other_compound:
        :param mappings:
        :return:
        """

        optimal_index = {1: -1, -1: -1, 2: -1, 0: -1}
        min_count =  {1: float("inf"), -1: float("inf"), 2: float("inf"), 0: float("inf")}
        for i, mapping in enumerate(mappings):
            one_stereo_counts, the_other_stereo_counts = self.compare_chemical_details_with_mapping(the_other_compound, mapping)
            if not one_stereo_counts and not the_other_stereo_counts:
                optimal_index[0] = i
                min_count[0] = 0
            if not one_stereo_counts:
                if one_stereo_counts < min_count[1]:
                    optimal_index[1] = i
                    min_count[1] = one_stereo_counts
            elif not the_other_stereo_counts:
                if the_other_stereo_counts < min_count[-1]:
                    optimal_index[-1] = i
                    min_count[-1] = the_other_stereo_counts
            else:
                if one_stereo_counts + the_other_stereo_counts < min_count[2]:
                    min_count[2] = one_stereo_counts + the_other_stereo_counts
                    optimal_index[2] = i
            # determine which one is the best.
        if min(min_count.values()) < float("inf"):
            relationship = self.determine_relationship(min_count)
            # return relationship and atom mappings
            final_mappings = mappings[optimal_index[relationship]]
            return relationship, { key: [final_mappings[key]] for key in final_mappings }
        return None, None

    @staticmethod
    def determine_relationship(unmapped_count):
        """
        To determine the relationship between two compounds when there are multiple possible atom mappings.
        We try to map as many details as possible. In other words, try to minimize the unmapped details.

        :param unmapped_count: the dictionary of relationship to the count of unmapped details.
        :type unmapped_count: :py:class:`dict`.
        :return: the relationship between the two compounds.
        :rtype: :py:class:`int`.
        """
        if unmapped_count[0] == 0:
            # equivalent mapping.
            return 0
        if unmapped_count[-1] < unmapped_count[1]:
            return -1
        elif unmapped_count[-1] > unmapped_count[1]:
            return -1
        elif unmapped_count[-1] == unmapped_count[1] and unmapped_count[-1] != float("inf"):
            return 1
        return 2

    def circular_pair_relationship(self, the_other_compound):
        """
        To determine the relationship of two compounds with interchangeable circular and linear representations.
        We first find the critical atoms that involve in the formation of ring. There can be several possibilities.
        Then we break the ring, and restore the double bond in the aldehyde group that forms the ring.
        Finally, check if the updated structure is the same with the other compound. And determine the relationship
        between the two compounds as well as generate the atom mappings.

        :param the_other_compound: the other :class:`~MDH.compound.Compound` entity.
        :type the_other_compound: :class:`~MDH.compound.Compound`.
        :return: the relationship and the atom mappings between the two compounds.
        :rtype: :py:class:`int` and :py:class:`dict`.
        """
        # default one compound should have a cycle.
        optimal_mappings = {1: None, -1: None, 2: None, 0: None}
        min_count = {1: float("inf"), -1: float("inf"), 2: float("inf"), 0: float("inf")}
        critical_atom_list = self.find_critical_atom_in_cycle()
        the_other_color = the_other_compound.backbone_color_identifier(r_groups=True) + \
                          the_other_compound.metal_color_identifier(details=False)
        for critical_atoms in critical_atom_list:
            self.break_cycle(critical_atoms)
            this_color = self.backbone_color_identifier(r_groups=True) + self.metal_color_identifier(details=False)
            if this_color == the_other_color:
                excluded_atoms_the_other = the_other_compound.exclude_atoms([self.atoms[i].color for i in critical_atoms])
                one_chemical_details = self.get_chemical_details(critical_atoms)
                the_other_chemical_details = the_other_compound.get_chemical_details(excluded_atoms_the_other)
                relationship, mis_count = self.compare_chemical_details(one_chemical_details, the_other_chemical_details)
                if mis_count < min_count[relationship]:
                    min_count[relationship] = mis_count
                    optimal_mappings[relationship] = self.generate_atom_mapping_by_atom_color(the_other_compound)
            self.restore_cycle(critical_atoms)
        if min(min_count.values()) < float("inf"):
            relationship = self.determine_relationship(min_count)
            return relationship, optimal_mappings[relationship]
        else:
            return None, None

    def exclude_atoms(self, colors):
        """

        :param colors:
        :return:
        """

        excluded_index = []
        for i, atom in enumerate(self.atoms):
            if atom.color in colors:
                excluded_index.append(i)
        return excluded_index

    def break_cycle(self, critical_atoms):
        """

        :param critical_atoms:
        :return:
        """
        atom_o, atom_c, atom_oo = critical_atoms
        # break the cycle
        self.atoms[atom_o].remove_neighbors([atom_c])
        self.atoms[atom_c].remove_neighbors([atom_o])
        # update the single bond to double bond,
        self.bond_lookup[(atom_c, atom_oo)].update_bond_type("2")
        self.color_compound(r_groups=True, bond_stereo=False, atom_stereo=False, resonance=False, isotope_resolved=False,
                            charge=False)

    def restore_cycle(self, critical_atoms):
        """

        :param critical_atoms:
        :return:
        """
        atom_o, atom_c, atom_oo = critical_atoms
        self.atoms[atom_o].add_neighbors([atom_c])
        self.atoms[atom_c].add_neighbors([atom_o])
        self.bond_lookup[(atom_c, atom_oo)].update_bond_type("1")
        self.color_compound(r_groups=True, bond_stereo=False, atom_stereo=False, resonance=False,
                            isotope_resolved=False,
                            charge=False)

    def find_critical_atom_in_cycle(self):
        """

        :return:
        """

        critical_atoms =[]
        for atom in self.atoms:
            if atom.in_cycle and atom.default_symbol == "O":
                for neighbor_index in atom.neighbors:
                    neighbor = self.atoms[neighbor_index]
                    for next_neighbor_index in neighbor.neighbors:
                        if next_neighbor_index != atom.atom_number and self.atoms[next_neighbor_index].default_symbol == "O" \
                                and self.bond_lookup[(next_neighbor_index, neighbor_index)].bond_type == "1":
                            critical_atoms.append([atom.atom_number, neighbor_index, next_neighbor_index])
        return critical_atoms

    def update_atom_symbol(self, index, updated_symbol):
        """

        :param index:
        :param updated_symbol:
        :return:
        """

        for i in index:
            self.atoms[i].update_symbol(updated_symbol)

    def validate_mapping_with_r(self, the_other_compound, one_rs, mapping):
        """

        :param the_other_compound:
        :param one_rs:
        :param mapping:
        :return:
        """
        # this is to compare the R GROUPS.
        # self is the more generic one, which suggests it should have more linkages.
        # self is the subset.
        reverse_index = { mapping[key]: key for key in mapping }
        one_r_linkages = collections.Counter()
        for idx in one_rs:
            r_atom = self.atoms[idx]
            for neighbor_index in r_atom.neighbors:
                bond = self.bond_lookup[(idx, neighbor_index)]
                # find the R group bonded atom.
                one_r_linkages["{0}-{1}".format(reverse_index[neighbor_index], bond.bond_type)] += 1

        the_other_r_linkages = collections.Counter()
        for idx, atom in enumerate(the_other_compound.atoms):
            if idx in mapping.values():
                for neighbor_index in atom.neighbors:
                    if the_other_compound.atoms[neighbor_index].default_symbol != "H" and neighbor_index not in mapping.values():
                        bond = the_other_compound.bond_lookup[(idx, neighbor_index)]
                        the_other_r_linkages["{0}-{1}".format(idx, bond.bond_type)] += 1

        if all(one_r_linkages[x] >= the_other_r_linkages[x] for x in the_other_r_linkages):
            return True
        return False

    def compare_chemical_details_with_mapping(self, the_other_compound, mapping):
        """

        :param the_other_compound:
        :param mapping:
        :return:
        """
        one_stereo_counts, the_other_stereo_counts = 0, 0
        one_consistent_atoms, the_other_consistent_atoms = set(), set()
        # check the atom stereochemistry
        for one_atom_index in mapping:
            one_atom = self.atoms[one_atom_index]
            the_other_atom = the_other_compound.atoms[mapping[one_atom_index]]
            if one_atom.color_layers[1] == the_other_atom.color_layers[1]:
                one_consistent_atoms.add(one_atom.atom_number)
                the_other_consistent_atoms.add(the_other_atom.atom_number)
                if one_atom.atom_stereo_parity == the_other_atom.atom_stereo_parity:
                    continue
                if one_atom.atom_stereo_parity != "0":
                    one_stereo_counts += 1
                if the_other_atom.atom_stereo_parity != "0":
                    the_other_stereo_counts += 1
        # check bond stereochemistry
        for bond in self.bonds:
            atom_1, atom_2 = bond.first_atom_number, bond.second_atom_number
            if atom_1 in one_consistent_atoms and atom_2 in one_consistent_atoms:
                the_other_atom_1, the_other_atom_2 = mapping[atom_1], mapping[atom_2]
                the_other_bond = the_other_compound.bond_lookup[(the_other_atom_1, the_other_atom_2)]
                if bond.bond_stereo == the_other_bond.bond_stereo:
                    continue
                if bond.bond_stereo != "0":
                    one_stereo_counts += 1
                if the_other_bond.bond_stereo != "0":
                    the_other_stereo_counts += 1
        return one_stereo_counts, the_other_stereo_counts

    def optimal_mapping_with_r(self, the_other_compound, one_rs, mappings):
        """

        :param the_other_compound:
        :param one_rs:
        :param mappings:
        :return:
        """

        optimal_index = {1: -1, -1: -1, 2: -1, 0: -1}
        min_count = {1: float("inf"), -1: float("inf"), 2: float("inf"), 0: float("inf")}
        for i, mm in enumerate(mappings):
            if self.validate_mapping_with_r(the_other_compound, one_rs, mm):
                one_stereo_counts, the_other_stereo_counts = self.compare_chemical_details_with_mapping(
                    the_other_compound, mm)
                if one_stereo_counts == 0 and the_other_stereo_counts < min_count[-1]:
                    min_count[-1] = the_other_stereo_counts
                    optimal_index[-1] = i
                elif one_stereo_counts != 0 and one_stereo_counts + the_other_stereo_counts < min_count[2]:
                    min_count[2] = one_stereo_counts + the_other_stereo_counts
                    optimal_index[2] = i
        if min(min_count.values()) < float("inf"):
            relationship = self.determine_relationship(min_count)
            return relationship, mappings[optimal_index[relationship]]
        else:
            return None, None

    def with_r_pair_relationship(self, the_other_compound):
        """

        :param the_other_compound:
        :return:
        """
        # self is the substructure, more generic, contain less chemical details.
        one_rs = list(self.r_groups)
        the_other_rs = list(the_other_compound.r_groups)
        self.update_atom_symbol(one_rs, "H")
        the_other_compound.update_atom_symbol(the_other_rs, "H")
        self.color_compound(r_groups=True, atom_stereo=False, bond_stereo=False)
        the_other_compound.color(r_groups=True, atom_stereo=False, bond_stereo=False)
        one_to_one_mappings = self.find_mappings(the_other_compound, resonance=False, r_distance=True)
        # here we need to consider the r_distance atom color identifier, so we need to color compounds.
        relationship, optimal_mappings = self.optimal_mapping_with_r(the_other_compound, one_rs, one_to_one_mappings)
        if optimal_mappings:
            self.update_atom_symbol(one_rs, "R")
            the_other_compound.update_atom_symbol(the_other_rs, "R")
            return relationship, self.map_r_correspondents(one_rs, the_other_compound, optimal_mappings)

        # match loosely.
        self.color_compound(r_groups=True, bond_stereo=False, atom_stereo=False, resonance=True)
        the_other_compound.color_compound(r_groups=True, bond_stereo=False, atom_stereo=False, resonance=True)
        one_to_one_mappings = self.find_mappings(the_other_compound, resonance=True, r_distance=True)
        relationship, optimal_mappings = self.optimal_mapping_with_r(the_other_compound, one_rs, one_to_one_mappings)
        self.update_atom_symbol(one_rs, "R")
        the_other_compound.update_atom_symbol(the_other_rs, "R")
        if optimal_mappings:
            return relationship, self.map_r_correspondents(one_rs, the_other_compound, optimal_mappings)
        return None, None

    def map_r_correspondents(self, one_rs, the_other_compound, mappings):
        # again, the self compound is substructure, we need to figure out the R group atom in self compound and its
        # corresponding atoms in the other compound.
        """

        :param one_rs:
        :param the_other_compound:
        :param mappings:
        :return:
        """

        full_mappings = collections.defaultdict(list)
        for idx in one_rs:
            r_atom = self.atoms[idx]
            visited = set(mappings[neighbor_index] for neighbor_index in r_atom.neighbors)
            r_correspondents = [neighbor_index for index in visited for neighbor_index in the_other_compound.atoms[index].neighbors
                                if neighbor_index not in mappings.values()]
            visited |= set(r_correspondents)
            while r_correspondents:
                full_mappings[idx].extend(r_correspondents)
                r_correspondents = the_other_compound.get_next_layer_neighbors(r_correspondents, visited)
        for idx in mappings:
            full_mappings[idx].append(mappings[idx])
        return full_mappings

