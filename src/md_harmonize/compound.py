#!/usr/bin/python3

"""
md_harmonize.compound
~~~~~~~~~~~~~~~~~~~~~

This module provides the :class:`~md_harmonize.compound.Atom` class, the :class:`~md_harmonize.compound.Bond` class,
and the :class:`~md_harmonize.compound.Compound` class to construct a compound entity. Most of the instance
variables of these three classes are based on CTFile fields.

"""

import collections
import numpy
import heapq
from typing import *
from pathlib import Path
import ctfile
import itertools 
from .supplement import not_r_groups
from .supplement import standard_bond_counts
from .supplement import atomic_weights
from .supplement import metal_symbols
from .supplement import index_to_charge
from . import tools
from datetime import datetime
try:
    from . import BASS
except ImportError:
    from . import uncythonized_BASS as BASS


class Atom:

    """ Atom class describes the :class:`~md_harmonize.compound.Atom` entity in the compound. """

    def __init__(self, atom_symbol: str, atom_number: int, x: str = "0", y: str = "0", z: str = "0",
                 mass_difference: str = "0", charge: str = "0", atom_stereo_parity: str = "0",
                 hydrogen_count: str = "0", stereo_care_box: str = "0", valence: str = "0", h0designator: str = "0",
                 atom_atom_mapping_number: str = "0", inversion_retention_flag: str = "0", exact_change_flag: str = "0",
                 kat: str = "", in_cycle: bool = False):
        """Atom initializer.

        :param atom_symbol: atom_symbol.
        :param atom_number: atom_number.
        :param x: the atom x coordinate.
        :param y: the atom y coordinate.
        :param z: the atom z coordinate.
        :param mass_difference: difference from mass in periodic table.
        :param charge: charge.
        :param atom_stereo_parity: atom stereo parity.
        :param hydrogen_count: hydrogen_count.
        :param stereo_care_box: stereo_care_box.
        :param valence: valence.
        :param h0designator: h0designator (obsolete CTFile parameter).
        :param atom_atom_mapping_number: atom_atom_mapping_number.
        :param inversion_retention_flag: inversion_retention_flag.
        :param exact_change_flag: exact_change_flag.
        :param kat: KEGG atom type.
        :param in_cycle: whether the atom is in cycle.
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
        self.in_cycle = in_cycle
        self.bond_counts = 0
        self.group_id = 0
        self.double_bond_counts = 0
        self.distance_to_r = 0
        self.kat = kat
        self.is_r = True if ("A" in self.atom_symbol or "R" in self.atom_symbol or "*" in self.atom_symbol or
                             self.atom_symbol == "X") and self.atom_symbol not in not_r_groups else False
        self.default_symbol = "R" if self.is_r else self.atom_symbol

    def update_symbol(self, symbol: str) -> str:
        """
        To update the atom symbol.

        :param symbol: the updated atom symbol.
        :return: the updated atom_symbol.
    """
        self.atom_symbol = self.default_symbol = symbol
        return self.atom_symbol

    def update_atom_number(self, index: int) -> int:
        """
        To update the atom number.

        :param index: the updated atom number.
        :return: the updated atom number.
    """
        self.atom_number = index
        return self.atom_number

    def remove_neighbors(self, neighbors: list) -> list:
        """
        To remove neighbors from the atom.

        :param neighbors: the list of neighbors that will be removed from this atom.
        :return: the updated list of neighbors.
    """
        for atom_index in neighbors:
            if atom_index in self.neighbors:
                self.neighbors.remove(atom_index)
        return self.neighbors

    def add_neighbors(self, neighbors: list) -> list:
        """
        To add neighbors to the atom.

        :param neighbors: the list of neighbors that will be added to this atom.
        :return: the updated list of neighbors.
    """
        for atom_index in neighbors:
            if atom_index not in self.neighbors:
                self.neighbors.append(atom_index)
        return self.neighbors

    def update_stereochemistry(self, stereo: str) -> str:
        """
        To update the atom stereochemistry.

        :param stereo: the updated atom stereochemistry.
        :return: the updated atom stereochemistry.
    """
        self.atom_stereo_parity = stereo
        return self.atom_stereo_parity
    
    def color_atom(self, isotope_resolved: bool = False, charge: bool = False, atom_stereo: bool = False) -> str:
        """
        To generate the atom color of the zero layer.

        :param isotope_resolved: If true, add isotope information when constructing colors.
        :param charge: If true, add charge information when constructing colors.
        :param atom_stereo: If true, add atom stereochemistry information when constructing colors.
        :return: the atom color of the zero layer.
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

    def reset_color(self) -> None:
        """
        Reset the atom color.

        :return: None.
    """
        self.color_0 = ""
        self.color = ""
        self.color_layers = {}
    
    def update_kat(self, kat: str) -> str:
        """
        To update the atom KEGG atom type.

        :param kat: the KEGG atom type for this atom,
        :return: the updated KEGG atom type.
    """
        self.kat = kat
        return self.kat

    def update_cycle(self, cycle_status: bool) -> bool:
        """
        To update the cycle status of the atom

        :param cycle_status: whether the atom is in cycle
        :return: cycle status
        """
        self.in_cycle = cycle_status
        return self.in_cycle

    def clone(self):
        """
        To clone the atom.

        :return: the cloned atom.
        """
        return Atom(self.atom_symbol, self.atom_number, x=str(self.x), y=str(self.y), z=str(self.z),
                    mass_difference=self.mass_difference, charge=str(self.charge),
                    atom_stereo_parity=self.atom_stereo_parity, hydrogen_count=self.hydrogen_count,
                    stereo_care_box=self.stereo_care_box, valence=self.valence, h0designator=self.h0designator,
                    atom_atom_mapping_number=self.atom_atom_mapping_number,
                    inversion_retention_flag=self.inversion_retention_flag, exact_change_flag=self.exact_change_flag,
                    kat=self.kat, in_cycle=self.in_cycle)


class Bond:

    """ Bond class describes the :class:`~md_harmonize.compound.Bond` entity in the compound. """

    def __init__(self, first_atom_number: str, second_atom_number: str, bond_type: str, bond_stereo: str = "0",
                 bond_topology: str = "0", reacting_center_status: str = "0"):
        """Bond initializer.

        :param first_atom_number: the index of the first atom forming this bond.
        :param second_atom_number: the index of the second atom forming this bond.
        :param bond_type: the bond type. (1 = Single, 2 = Double, 3 = Triple, 4 = Aromatic, 5 = Single or Double, 6 = Single or Aromatic, 7 = double or Aromatic 8 = Any)
        :param bond_stereo: the bond stereo. (Single bonds: 0 = not stereo, 1 = Up, 4 = Either, 6 = Down; Double bonds: determined by x, y, z coordinates)
        :param bond_topology: bond topology. (O = Either, 1 = Ring, 2 = Chain)
        :param reacting_center_status: reacting center status.
    """

        self.first_atom_number = int(first_atom_number) - 1
        self.second_atom_number = int(second_atom_number) - 1
        self.bond_type = bond_type.strip()
        self.bond_stereo = bond_stereo.strip() if bond_stereo.strip() != "4" and bond_stereo.strip() != "8" else "0"
        self.bond_topology = bond_topology.strip()
        self.reacting_center_status = reacting_center_status.strip()
 
    def update_bond_type(self, bond_type: str) -> str:
        """
        To update the bond type.

        :param bond_type: the updated bond type.
        :return: the updated bond type.
    """
        self.bond_type = str(bond_type)
        return self.bond_type

    def update_stereochemistry(self, stereo: str) -> str:
        """
        To update the bond stereochemistry.

        :param stereo: the updated bond stereochemistry.
        :return: the updated bond stereochemistry.
    """
        self.bond_stereo = str(stereo)
        return self.bond_stereo

    def update_first_atom(self, index: int) -> int:
        """
        To update the first atom number of the bond.

        :param index: the updated first atom number.
        :return: the updated first atom number.
    """
        self.first_atom_number = index
        return self.first_atom_number

    def update_second_atom(self, index: int) -> int:
        """
        To update the second atom number of the bond.

        :param index: the updated second atom number.
        :return: the updated second atom number.
    """
        self.second_atom_number = index
        return self.second_atom_number

    def clone(self):
        """
        To clone the bond.

        :return: the cloned bond.
    """
        return Bond(str(self.first_atom_number+1), str(self.second_atom_number+1), self.bond_type,
                    bond_stereo=self.bond_stereo, bond_topology=self.bond_topology,
                    reacting_center_status=self.reacting_center_status)


class Compound:

    """ Compound class describes the :class:`~md_harmonize.compound.Compound` entity. """

    def __init__(self, compound_name: str, atoms: list, bonds: list) -> None:
        """Compound initializer.

        :param compound_name: the compound name.
        :param atoms: a list of :class:`~md_harmonize.compound.Atom` entities in the compound.
        :param bonds: a list of :class:`~md_harmonize.compound.Bond` entities in the compound.
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
        self._distance_matrix = None

    def __str__(self):
        """
        To generate the string representation of the compound.

        """
        str = ""
        for atom in self.atoms:
            str += "{0} {1}\n".format(atom.atom_number, atom.default_symbol)
        for bond in self.bonds:
            str += "{0} {1} {2}\n".format(bond.first_atom_number, bond.second_atom_number, bond.bond_type)
        return str

    def encode(self) -> tuple:
        """
        To clone the compound.

        :return: the cloned compound.
    """
        return self.name, [atom.clone() for atom in self.atoms], [bond.clone() for bond in self.bonds]

    @property
    def name(self) -> str:
        """
        To get the compound name.

        :return: the compound name.
    """
        return self.compound_name

    @staticmethod
    def molfile_name(molfile: str):
        """
        Create the compound entity based on the molfile representation.

        :param molfile: the filename of the molfile.
        :return: the constructed compound entity.
    """
        compound_name = Path(molfile).stem
        try:
            with open(molfile, 'r') as infile:
                ct_object = ctfile.load(infile)
            atoms = [Atom(atom.atom_symbol, i, atom['x'], atom['y'], atom['z'],  atom['mass_difference'], atom.charge,
                           atom['atom_stereo_parity'], atom['hydrogen_count'], atom['stereo_care_box'], atom['valence'],
                           atom['h0designator'], atom['atom_atom_mapping_number'], atom['inversion_retention_flag'],
                           atom['exact_change_flag']) for i, atom in enumerate(ct_object.atoms) ]
            bonds = [Bond(bond['first_atom_number'], bond['second_atom_number'], bond['bond_type'], bond['bond_stereo'],
                          bond['bond_topology'], bond['reacting_center_status']) for bond in ct_object.bonds]
            return Compound(compound_name, atoms, bonds)
        except Exception as e:
            print(e)
            print(compound_name + " cannot not be converted to ctifle object")
            return None

    @property
    def formula(self) -> str:
        """
        To construct the formula of this compound (only consider heavy atoms).

        :return: string formula of the compound.
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
    def composition(self) -> dict:
        """
        To get the atom symbols and bond types in the compound.

        :return: the atom and bond information of the compound
    """
        atom_composition = collections.Counter([atom.default_symbol for atom in self.atoms])
        bond_composition = collections.Counter([bond.bond_type for bond in self.bonds])
        atom_composition.update(bond_composition)
        return atom_composition

    @property
    def r_groups(self) -> list:
        """
        To get all the R groups in the compound.

        :return: the list of index of all the R groups.
    """
        return [index for index, atom in enumerate(self.atoms) if atom.is_r]

    def contains_r_groups(self) -> bool:
        """
        To check if the compound contains R group(s).

        :return: bool whether the compound contains R group.
    """
        return self.r_groups != []

    def has_isolated_atoms(self) -> bool:
        """
        To check if the compound has atoms that have no connections to other atoms.

        :return: bool whether the compound has isolated atoms.
    """
        for atom in self.atoms:
            if atom.bond_counts == 0:
                return True
        return False

    @property
    def metal_index(self) -> list:
        """
        To get the metal elements in the compound.

        :return: a list of atom numbers of metal elements.
    """
        return [index for index, atom in enumerate(self.atoms) if atom.default_symbol in metal_symbols]
    
    @property
    def h_index(self) -> list:
        """
        To get all H in the compound.

        :return: a list of atom numbers corresponding to H.
    """
        return [index for index, atom in enumerate(self.atoms) if atom.default_symbol == "H"]

    @property
    def heavy_atoms(self) -> list:
        """
        To get all the heavy atoms in the compound.

        :return: a list of atom numbers corresponding to heavy atoms.
    """
        return [atom for atom in self.atoms if atom.default_symbol != "H"]

    @property
    def index_of_heavy_atoms(self) -> dict:
        """
        To map the atom number to index in the heavy atom list.

        :return: the dictionary of atom number to atom index of heavy atoms.
    """
        return { atom.atom_number: i for i, atom in enumerate(self.heavy_atoms) }

    # def is_symmetric(self, atom_index):
    #
    #     for index in atom_index:
    #         color = self.atoms[index].color
    #         if len(self.color_compound()[color]) > 1:
    #             return False
    #     return True
       
    def color_groups(self, excluded=None) -> dict:
        """
        To update the compound color groups after coloring.

        :return: the dictionary of atom color with the list of atom number.
    """
        if not excluded:
            excluded = []
        color_groups = collections.defaultdict(list)
        for atom in self.atoms:
            if atom.atom_number not in excluded:
                color_groups[atom.color].append(atom.atom_number)
        return color_groups

    def detect_abnormal_atom(self) -> dict:
        """
        To find the atoms with invalid bond counts.

        :return: a list of atom numbers with invalid bond counts.
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

    def curate_invalid_n(self) -> None:
        """
        To curate the charge of invalid N atoms.

        :return: None.
    """
        abnormal_atoms = self.detect_abnormal_atom()
        for atom_index in abnormal_atoms["N"]:
            atom = self.atoms[atom_index]
            atom.charge += 1
            for neighbor in atom.neighbors:
                if self.atoms[neighbor].default_symbol == "O" and \
                        self.bond_lookup[(atom.atom_number, neighbor)].bond_type == "2":
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

    def update_aromatic_bond_type(self, cycles: list) -> None:
        """
        Update the aromatic bond types.
        Two cases: 1) change the bond in the aromatic ring to aromatic bond (bond type = 4); 2) change the double bond connecting to the aromatic ring to single bond.

        :param cycles: the list of cycles represented by aromatic atom index.
        :return: None.
    """
        atom_in_cycle = [atom for cycle in cycles for atom in cycle]
        for cycle in cycles:
            aromatic_bonds = self.extract_aromatic_bonds(cycle)
            for bond in aromatic_bonds:
                bond.update_bond_type("4")
        bond_out_of_cycle = self.extract_double_bond_connecting_cycle(atom_in_cycle)
        for bond in bond_out_of_cycle:
            bond.update_bond_type("1")

    # def update_aromatic_bond(self, aromatic_bonds, aromatic_atoms):
    #    """
    #    To update the bond type of aromatic bonds in the compound.
    #    Two steps are involved:
    #        1) Update the bond in the aromatic ring to type "4".
    #        2) Update the double bond connecting to the aromatic ring to "1".
    #
    #    :param aromatic_bonds: a list of aromatic bonds in the compound (represented by the first and second atom number).
    #    :type aromatic_bonds: :py:class:`list`.
    #    :param aromatic_atoms: a list of aromatic atoms in the aromatic ring.
    #    :type aromatic_atoms: :py:class:`list`.
    #    :return: None.
    #    :rtype: :py:obj:`None`.
    #    """
    #
    #    for atom_i, atom_j in aromatic_bonds:
    #        self.bond_lookup[(atom_i, atom_j)].update_bond_type("4")
    #    bond_out_of_cycle = self.extract_double_bond_connecting_cycle(aromatic_atoms)
    #    for bond in bond_out_of_cycle:
    #        bond.update_bond_type("1")

    def extract_double_bond_connecting_cycle(self, atom_in_cycle: list) -> list:
        """
        To extract the double bonds connecting to the atom in the aromatic cycles.

        :param atom_in_cycle: the list of aromatic cycles represented by aromatic atom index.
        :return: the list of outside double bond connecting to the atom in the aromatic cycles.
    """
        double_bond_connecting_cycle = []
        for atom_index in atom_in_cycle:
            if self.atoms[atom_index].default_symbol == "C":
                for neighbor_index in self.atoms[atom_index].neighbors:
                    if neighbor_index not in atom_in_cycle and \
                            self.bond_lookup[(atom_index, neighbor_index)].bond_type == "2":
                        double_bond_connecting_cycle.append(self.bond_lookup[(atom_index, neighbor_index)])
        return double_bond_connecting_cycle

    def extract_aromatic_bonds(self, cycle: list) -> list:
        """
        Extract the aromatic bonds based on the atoms in the cycle.

        :param cycle: the list of aromatic cycles represented by aromatic atom index.
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

    def separate_connected_components(self, bonds: Union[list, set]) -> list:
        """
        This is used in constructing the aromatic substructures detected by the Indigo method.
        A compound can have several disjoint aromatic substructures. Here, we need to find the disjoint parts.
        The basic idea is union-find. We union atoms that are connected by a bond.

        :param bonds: the list of bonds representing by the atom numbers forming the bond.
        :return: a list of separate components represented by a list atom numbers in the component.
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

    def connected_components(self) -> dict:
        """
        Detect the connected components in the compound structure (using the breadth first search).
        Cases when not all the atoms are connected together.

        :return: the dictionary of the connected components.
    """
        visited = set()
        group_id = 0
        components = collections.defaultdict(list)
        for index, atom in enumerate(self.atoms):
            if index not in visited:
                ques = collections.deque([index])
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

    def calculate_distance_to_r_groups(self) -> None:
        """
        To calculate the distance of each atom to its nearest R group (using the dijkstra's algorithm).

        :return: None:
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
                    distance_matrix[self.index_of_heavy_atoms[atom_index]] = \
                        min(distance_matrix[self.index_of_heavy_atoms[atom_index]], dist)
                    for neighbor in self.atoms[atom_index].neighbors:
                        if self.atoms[neighbor].default_symbol != "H":
                            if neighbor not in seen or seen[neighbor] > dist+ 1:
                                if distance_matrix[self.index_of_heavy_atoms[neighbor]] > dist+1:
                                    heapq.heappush(ques, (dist+1, neighbor))
        
            for i, dist in enumerate(distance_matrix):
                self.heavy_atoms[i].distance_to_r = dist

    def find_cycles(self, short_circuit: bool = False, cutoff: int = 40, seconds=50) -> list:
        """
        To find the cycles in the compound.

        :param short_circuit: whether to take short path.
        :param cutoff: limit of cycle length.
        :param seconds: the timeout limit.
        :return: the list of cycles in the compound.
        """
        try:
            with tools.timeout_context(seconds=seconds):
                return self.find_cycles_helper(short_circuit=short_circuit, cutoff=cutoff)
        except Exception as e:
            print("Cycles in compound {0} can hardly be detected: {1}".format(self.name, e))
            pass
        return []

    # def find_cycles(self, short_circuit: bool = False, cutoff: int = 40, seconds: int = 50) -> None:
    #     """
    #     To find the cycles in the compound with timeout limit.
    #
    #     :param short_circuit: whether to take short path.
    #     :param cutoff: limit of cycle length.
    #     :param seconds: the timeout limit.
    #     """
    #     try:
    #         tools.timeout(self.find_cycles_helper, (short_circuit, cutoff,), seconds=seconds)
    #         print("in timeout find_cycles", self.has_cycle)
    #         print("in timeout find_cycles", [atom.atom_number for atom in self.atoms if atom.in_cycle])
    #     except Exception as e:
    #         print("Cycles in compound {0} can hardly be detected: {1}".format(self.name, e))
    #         pass

    def find_cycles_helper(self, short_circuit: bool = False, cutoff: int = 40) -> list:
        """
        Executing function to find the cycles in the compound.

        :param short_circuit: whether to take short path.
        :param cutoff: limit of cycle length.
        :return: the list of cycles in the compound
        """
        not_cyclic, cyclic, all_cycles = [], [], []

        def prune() -> None:
            while 1:
                terminate = True
                for prune_atom in self.atoms:
                    if prune_atom.atom_number not in not_cyclic:
                        # for atom with cycle, it at least has 2 neighbors to form a cycle.
                        if len([neighbor for neighbor in prune_atom.neighbors if neighbor not in not_cyclic]) < 2:
                            not_cyclic.append(prune_atom.atom_number)
                            terminate = False
                if terminate:
                    break

        def search(start_index: int) -> list:

            # the paths is a list of list, which means each path is a list.
            paths, cycles = collections.deque([[start_index]]), []
            while paths:
                indices = paths.popleft()
                for neighbor in self.atoms[indices[-1]].neighbors:
                    if neighbor not in not_cyclic:
                        if neighbor not in indices and len(indices) < cutoff:
                            paths.append(list(indices) + [self.atoms[neighbor].atom_number])
                        elif neighbor == start_index and len(indices) != 2:
                            cycles.append(list(indices) + [self.atoms[neighbor].atom_number])
                            if short_circuit:
                                return cycles
            if cycles:
                return cycles
        prune()
        for atom in self.atoms:
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
                    self.atoms[index].update_cycle(True)
        if all_cycles:
            self.has_cycle = True
        return [list(x) for x in set(tuple(x) for x in [sorted(i[:-1]) for l in all_cycles for i in l])]

    def structure_matrix(self, resonance: bool = False, backbone: bool = False) -> numpy.ndarray:
        """
        To construct graph structural matrix of this compound.
        matrix[i][j] = 0 suggests the two atoms are not connected directly.
        Other integer represented the bond type connecting the two atoms.

        :param resonance: bool whether to ignore the difference between single and double bonds.
        :param backbone: bool whether to ignore bond types. This is for parsing atoms mappings from KEGG RCLASS.
        :return: the constructed structure matrix for this compound.
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
    def distance_matrix(self) -> numpy.ndarray:
        """
        To construct the distance matrix of the compound (using the Floyd Warshall Algorithm).
        distance[i][j] suggests the distance between atom i and j.

        :return: the distance matrix of the compound.
    """
        if self._distance_matrix is None:
            if self.heavy_atoms:
                distance_matrix = numpy.ones((len(self.heavy_atoms), len(self.heavy_atoms)),dtype=numpy.uint16 )\
                                  * len(self.heavy_atoms)
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
                self._distance_matrix = distance_matrix
        return self._distance_matrix
    
    def update_color_tuple(self, resonance: bool = False) -> None:
        """
        To update the color tuple of the atoms in the compound. This color tuple includes information of its neighboring
        atoms and bonds. Here, we don't need to consider backbone since this part was initially designed for aromatic
        substructure detection and only double and single bonds are considered.

        :param resonance: bool whether to ignore the difference between single and double bonds.
        :return: None.
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
    
    def find_mappings(self, the_other, resonance: bool = True, r_distance: bool = False,
                      backbone: bool = False) -> list:
        """
        Find the one to one atom mappings between two compounds using the BASS algorithm. The other compound  is supposed be
        contained in the self compound.

        :param the_other: the mappings compound entity.
        :param resonance: whether to ignore the difference between single and double bonds.
        :param r_distance: whether to take account of the position of R groups.
        :param backbone: whether to ignore the bond types.
        :return: the list of atom mappings in the heavy atom order.
    """
        self.update_color_tuple(resonance=resonance)
        the_other.update_color_tuple(resonance=resonance)
        mappings = []
        mapping_matrix = BASS.make_mapping_matrix(the_other, self, True, True, r_distance)
        if mapping_matrix is not None:
            mappings = BASS.find_mappings(the_other.structure_matrix(resonance=resonance, backbone=backbone),
                                          the_other.distance_matrix, self.structure_matrix(resonance=resonance,
                                                                                           backbone=backbone),
                                          self.distance_matrix, mapping_matrix)
        # for the mappings, the from_idx, to_idx in enumerate(mapping), from_idx is in the_other_compound, to_idx is in
        # the self.
        one_to_one_mappings = []
        for sub in mappings:
            cur_mappings = {}
            for from_index, to_index in enumerate(sub):
                cur_mappings[self.heavy_atoms[to_index].atom_number] = the_other.heavy_atoms[from_index].atom_number
            one_to_one_mappings.append(cur_mappings)
        return one_to_one_mappings

    def find_mappings_reversed(self, the_other, resonance: bool = True, r_distance: bool = False,
                               backbone: bool = False) -> list:
        """
        Find the one to one atom mappings between two compounds using the BASS algorithm. The self compound is supposed to
        be contained in the other compound.

        :param the_other: the mappings compound entity.
        :param resonance: whether to ignore the difference between single and double bonds.
        :param r_distance: whether to take account of the position of R groups.
        :param backbone: whether to ignore the bond types.
        :return: the list of atom mappings in the heavy atom order.
    """
        self.update_color_tuple(resonance=resonance)
        the_other.update_color_tuple(resonance=resonance)
        mappings = []
        mapping_matrix = BASS.make_mapping_matrix(self, the_other, True, True, r_distance)
        if mapping_matrix is not None:
            mappings = BASS.find_mappings(self.structure_matrix(resonance=resonance, backbone=backbone),
                                          self.distance_matrix,
                                          the_other.structure_matrix(resonance=resonance, backbone=backbone),
                                          the_other.distance_matrix,
                                          mapping_matrix)
        # for the mappings, the from_idx, to_idx in enumerate(mapping), from_idx is in the_other_compound, to_idx is in
        # the self.
        one_to_one_mappings = []
        for sub in mappings:
            cur_mappings = {}
            for from_index, to_index in enumerate(sub):
                cur_mappings[self.heavy_atoms[from_index].atom_number] = the_other.heavy_atoms[to_index].atom_number
            one_to_one_mappings.append(cur_mappings)
        return one_to_one_mappings

    def map_resonance(self, the_other, r_distance: bool = False, seconds: int = 50) -> list:
        """
        Check if the resonant mappings are valid between two compound structures.

        :param the_other: the mappings compound entity.
        :param r_distance: to take account of the position of R groups.
        :param seconds: the timeout limit.
        :return: the list of valid atom mappings between the two compound structures.
        """
        try:
            mapping = tools.timeout(self.map_resonance_helper, (the_other, r_distance,), seconds=seconds)
            return mapping
        except Exception as exception:
            return []

    def map_resonance_helper(self, the_other, r_distance: bool = False) -> list:
        """
        Check if the resonant mappings are valid between the two compound structures. If the mapped atoms don't share
        the same local coloring identifier, we check if the difference is caused by the position of double bonds.
        Find the three atoms involved in the resonant structure and check if one of the atoms is not C.
        N (a)            N (a)
        / \\             // \
        (b) C   N (c)    (b) C   N (c)

        In addition, the self compound is supposed to be more generic, which means has fewer atoms. Therefore, atoms in
        self compound can all be mapped to the other compound.

        :param the_other: the mappings compound entity.
        :param r_distance: to take account of position of R groups.
        :return: the list of valid atom mappings between the two compound structures.
    """
        self.color_compound(r_groups=True, atom_stereo=False, bond_stereo=False)
        the_other.color_compound(r_groups=True, atom_stereo=False, bond_stereo=False)

        one_to_one_mappings = self.find_mappings(the_other, resonance=True, r_distance=r_distance)

        valid_mappings = []
        for cur_mappings in one_to_one_mappings:
            flag = True
            reversed_mappings = { cur_mappings[key] : key for key in cur_mappings }
            for from_index in sorted(cur_mappings.keys()):
                to_index = cur_mappings[from_index]

                # Detect if the two atoms have the same local binding environment by comparing the first layer color
                # identifier.
                if 1 in self.atoms[from_index].color_layers and 1 in the_other.atoms[to_index].color_layers and \
                        self.atoms[from_index].color_layers[1] == the_other.atoms[to_index].color_layers[1]:
                    continue
                three = {from_index}
                from_dbond_atom_index = self.find_double_bond_linked_atom(from_index)
                to_dbond_atom_index = the_other.find_double_bond_linked_atom(to_index)

                # cannot find directly linked double bond.
                if from_dbond_atom_index == -1 and to_dbond_atom_index == -1:
                    flag = False
                # Need to conduct second search
                # atom b in the picture
                elif from_dbond_atom_index == -1 and to_dbond_atom_index in reversed_mappings:
                    # here we can find the atom a
                    reversed_from_dbond_atom_index = reversed_mappings[to_dbond_atom_index]
                    three.add(reversed_from_dbond_atom_index)
                    from_next_dbond_atom_index = self.find_double_bond_linked_atom(reversed_from_dbond_atom_index) 
                    # based on atom a, we can find atom c
                    if from_next_dbond_atom_index != -1:
                        three.add(from_next_dbond_atom_index)

                # atom c in the picture
                elif to_dbond_atom_index == -1:
                    three.add(from_dbond_atom_index) # here we find atom a
                    reversed_to_dbond_atom_index = cur_mappings[from_dbond_atom_index]
                    # based on atom a find atom b in the other compound.
                    to_next_dbond_atom_index = the_other.find_double_bond_linked_atom(reversed_to_dbond_atom_index)
                    if to_next_dbond_atom_index != -1 and to_next_dbond_atom_index in reversed_mappings:
                        three.add(reversed_mappings[to_next_dbond_atom_index])

                # atom a in the picture
                # find the other two atoms directly.
                elif from_dbond_atom_index != -1 and to_dbond_atom_index != -1 and to_dbond_atom_index in reversed_mappings:
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

    def find_double_bond_linked_atom(self, i: int) -> int:
        """
        Find the atom that is doubly linked to the target atom i.

        :param i: the ith atom in the compound.
        :return: the index of the doubly linked atom.
    """
        for neighbor_index in self.atoms[i].neighbors:
            if neighbor_index in self.index_of_heavy_atoms and self.bond_lookup[(i, neighbor_index)].bond_type == "2":
                return neighbor_index
        return -1

    def define_bond_stereochemistry(self) -> None:
        """
        Define the stereochemistry of double bonds in the compound.

        :return: None.
    """
        for bond in self.bonds:
            if bond.bond_type == "2":
                bond.update_stereochemistry(self.calculate_bond_stereochemistry(bond))
    
    def calculate_bond_stereochemistry(self, bond: Bond) -> int:
        """
        Calculate the stereochemistry of the double bond based on its geometric properties. The line of the double bond divides
        the plane into two parts. For the atoms forming the double bond, it normally has two branches. If the two
        branches are not the same, we call them heavy side and light side (heavy side containing atoms with heavier
        atomic weights). We determine the bond stereochemistry by checking if the two heavy sides lie on the same part
        of the divided plane.
        
        :param bond: the bond entity.
        :return: the calculated bond stereochemistry.
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

        if (first_light_branch and first_light_branch == first_heavy_branch) or \
                (second_light_branch and second_light_branch == second_heavy_branch):
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
    def calculate_y_coordinate(slope: float, b: float, atom: Atom) -> float:
        """
        Calculate the y coordinate of the atom based on the linear function: y = slope * x + b

        :param slope: the slope of the targeted line.
        :param b: the intercept of the targeted line.
        :param atom: the atom entity.
        :return: the calculated y coordinate.
    """
        return atom.x * slope + b

    def collect_atomic_weights_of_neighbors(self, neighbors: list) -> list:
        """
        To collect the atomic weights of the current layer's neighbors.

        :param neighbors: the list of atom numbers of neighbors.
        :return: the list of atomic weights for this layer's neighbors.
    """
        neighbor_atomic_weights = [atomic_weights[self.atoms[index].default_symbol] for index in neighbors]
        neighbor_atomic_weights.sort(reverse=True)
        return neighbor_atomic_weights

    def compare_branch_weights(self, neighbors: list, atom_forming_double_bond: Atom) -> tuple:
        """
        To determine the heavy and light branches that connect to the atom forming the double bond. This is based on
        comparison of the atomic weights of the two branches (breadth first algorithm).

        :param neighbors: the list of atom numbers of the atoms that connect the atom forming the double bond.
        :param atom_forming_double_bond: the atom that forms the bond.
        :return: heavy and light branches. [heavy_side, light_side]
    """
        if len(neighbors) < 2:
            return neighbors[0], None

        one_atom, the_other_atom = neighbors[0], neighbors[1]

        if atomic_weights[self.atoms[one_atom].default_symbol] > atomic_weights[self.atoms[the_other_atom].default_symbol]:
            return one_atom, the_other_atom

        elif atomic_weights[self.atoms[one_atom].default_symbol] < atomic_weights[self.atoms[the_other_atom].default_symbol]:
            return the_other_atom, one_atom

        else:
            one_neighbors = [one_atom]
            one_visited = {atom_forming_double_bond.atom_number}
            the_other_neighbors = [the_other_atom]
            the_other_visited = {atom_forming_double_bond.atom_number}
            one_neighbor_atomic_weight_list = self.collect_atomic_weights_of_neighbors(one_neighbors)
            the_other_neighbor_atomic_weight_list = self.collect_atomic_weights_of_neighbors(the_other_neighbors)

            while one_neighbor_atomic_weight_list == the_other_neighbor_atomic_weight_list and \
                    one_neighbor_atomic_weight_list:
                one_neighbors = self.get_next_layer_neighbors(one_neighbors, one_visited)
                one_neighbor_atomic_weight_list = self.collect_atomic_weights_of_neighbors(one_neighbors)
                the_other_neighbors = self.get_next_layer_neighbors(the_other_neighbors, the_other_visited)
                the_other_neighbor_atomic_weight_list = self.collect_atomic_weights_of_neighbors(the_other_neighbors)

            if tuple(one_neighbor_atomic_weight_list) > tuple(the_other_neighbor_atomic_weight_list):
                return one_atom, the_other_atom
            elif tuple(one_neighbor_atomic_weight_list) < tuple(the_other_neighbor_atomic_weight_list):
                return the_other_atom, one_atom
            else:
                return one_atom, one_atom

    def get_next_layer_neighbors(self, cur_layer_neighbors: list, visited: set, excluded: list = None) -> list:
        """
        To get the next layer's neighbors.

        :param cur_layer_neighbors: the list of atom numbers of the current layer.
        :param visited: the atom numbers that have already been visited.
        :param excluded: the list of atom numbers that should not be included in the next layer.
        :return: the neighboring atom numbers of the next layer.
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

    def color_compound(self, r_groups: bool = True, bond_stereo: bool = False, atom_stereo: bool = False,
                       resonance: bool = False, isotope_resolved: bool = False, charge: bool = False,
                       backbone: bool = False) -> None:
        """
        To color the compound.

        :param r_groups:  If true, add R groups in the coloring.
        :param bond_stereo:  If true, add bond stereo detail when constructing colors.
        :param atom_stereo: If true, add atom stereo detail when constructing colors.
        :param resonance: If true, ignore the difference between double and single bonds.
        :param isotope_resolved: If true, add isotope detail when constructing colors.
        :param charge: If true, add charge detail when constructing colors.
        :param backbone: If true, ignore bond types in the coloring.
        :return: None.
    """
        excluded_index = self.metal_index + self.h_index
        excluded_index += self.r_groups if not r_groups else []
        atoms_to_color = [i for i in range(len(self.atoms)) if i not in excluded_index]
        self.reset_color()
        self.generate_atom_zero_layer_color(isotope_resolved=isotope_resolved, charge=charge, atom_stereo=atom_stereo)
        self.first_round_color(atoms_to_color, excluded_index=excluded_index, bond_stereo=bond_stereo, resonance=resonance,
                               backbone=backbone)
        self.curate_invalid_symmetric_atoms(atoms_to_color, excluded_index=excluded_index, bond_stereo=bond_stereo,
                                            resonance=resonance, backbone=backbone)
        self.color_metal(bond_stereo=bond_stereo, resonance=resonance)
    
    def reset_color(self) -> None:
        """
        To set the color of atoms in the compound to be empty.

        :return: None:
    """
        for atom in self.atoms:
            atom.reset_color()

    def generate_atom_zero_layer_color(self, isotope_resolved: bool = False, charge: bool = False,
                                       atom_stereo: bool = False) -> None:
        """
        To generate the color identifier of zero layer for each atom. We don't consider H and metals here.

        :param isotope_resolved: If true, add isotope detail when constructing colors.
        :param charge: If true, add charge detail when constructing colors.
        :param atom_stereo: If true, add atom stereochemistry detail when constructing colors.
        :return: None.
    """
        for index, atom in enumerate(self.atoms):
            atom.color_atom(isotope_resolved=isotope_resolved, charge=charge, atom_stereo=atom_stereo)

    def generate_atom_color_with_neighbors(self, atom_index: list, excluded: list = None, zero_core_color: bool = True,
                                           zero_neighbor_color: bool = True, resonance: bool = False,
                                           bond_stereo: bool = False, backbone: bool = False) -> dict:
        """
        To generate the atom color with its neighbors. We add this color name when we try to incorporate neighbors'
        information in naming.

        Here, we don't need to care about the atom stereo. It has been taken care of in generating color_0.

        Basic color formula: atom.color + [neighbor.color + bond.bond_type]

        :param atom_index: indices of atoms to color.
        :param excluded: the list of atom indices will be excluded from coloring.
        :param zero_core_color: If true, we use the atom.color_0 else atom.color for the core atom (first round coloring vs validation).
        :param zero_neighbor_color: If true, we use the atom.color_0 else atom.color for the neighbor atoms (first round coloring vs validation).
        :param resonance: If true, detect resonant compound pairs without distinguishing between double and single bonds.
        :param bond_stereo:  If true, add stereo detail of bonds when constructing colors.
        :param backbone: If true, ignore bond types in the coloring.
        :return: the dictionary of atom index and its color name.
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

    def first_round_color(self, atoms_to_color: list, excluded_index: list = None, bond_stereo: bool = False,
                          resonance: bool = False, backbone: bool = False, depth: int = 5000) -> None:
        """
        To do the first round of coloring this compound. We add neighbors' information layer by layer to the atom's
        color identifier until it has a unique identifier or all the atoms in the compound have been used for naming
        (based on the breadth first search algorithm).

        :param atoms_to_color: the list of atom numbers to be colored.
        :param excluded_index: the list of atom numbers to be excluded from coloring.
        :param bond_stereo: If true, add bond stereo detail when constructing colors.
        :param resonance: If true, ignore the difference between double and single bonds.
        :param backbone: If true, ignore bond types in the coloring.
        :param depth: the max depth of coloring.
        :return: None.
    """
        atom_color_with_neighbors = self.generate_atom_color_with_neighbors(atoms_to_color, excluded=excluded_index,
                                                                            zero_core_color=True,
                                                                            zero_neighbor_color=True,
                                                                            resonance=resonance,
                                                                            bond_stereo=bond_stereo,
                                                                            backbone=backbone)
        if not excluded_index:
            excluded_index = []

        atom_neighbors = collections.defaultdict(list)
        visited = collections.defaultdict(set)
        
        for index in atoms_to_color:
            atom_neighbors[index].append(index)
            visited[index].add(index)

        if depth == 5000:
            depth = len(atoms_to_color)

        i = 0
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
                    added += "({0}_{1})".format(color_elements[name], name) if color_elements[name] > 1 \
                        else "({0})".format(name)
                atom.color += added
                atom.color_layers[i+1] = atom.color
                atom_neighbors[atom.atom_number] = self.get_next_layer_neighbors(atom_neighbors[atom_index],
                                                                                 visited[atom_index],
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

    def invalid_symmetric_atoms(self, atoms_to_color: list, excluded_index: bool = None, bond_stereo: bool = False,
                                resonance: bool = False, backbone: bool = False) -> list:
        """
        To check if atoms with the same color identifier are symmetric.

        :param atoms_to_color: the list of atom numbers to be colored.
        :param excluded_index: the list of atom numbers to be excluded from coloring.
        :param bond_stereo: If true, add bond stereo detail when constructing colors.
        :param resonance: If true, ignore the difference between double and single bonds.
        :param backbone: If true, ignore bond types in the coloring.
        :return: the list of atom numbers to be recolored.
    """
        atom_color_with_neighbors = self.generate_atom_color_with_neighbors(atoms_to_color, excluded=excluded_index,
                                                                            zero_core_color=False,
                                                                            zero_neighbor_color=False,
                                                                            resonance=resonance,
                                                                            bond_stereo=bond_stereo,
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
                    target_color_list = [atom_color_with_neighbors[atom_index] for atom_index in
                                         current_layer[target_index] if atom_index not in excluded_index]
                    target_color_list.sort()
                    target_color = ''.join(target_color_list)

                    # check if neighbors of this layer are the same among all the atoms with the same color identifier.
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
                    # If they share the same neighbors, check the next layer.
                    if flag:
                        for atom_index in current_layer.keys():
                            next_layer_neighbors = self.get_next_layer_neighbors(current_layer[atom_index],
                                                                                 visited[atom_index],
                                                                                 excluded=excluded_index)
                            current_layer[atom_index] = next_layer_neighbors
                        # pruning check, to see if they have the same number of next layer's neighbors.
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

    def curate_invalid_symmetric_atoms(self, atoms_to_color: list, excluded_index: list = None,
                                       bond_stereo: bool = False, resonance: bool = False, backbone: bool = False) -> None:
        """
        To curate the color identifiers of invalid symmetric atoms.
        We recolor those invalid atoms using the full color identifiers of its neighbors layer by layer until the
        difference can be captured.

        :param atoms_to_color: the list of atom numbers of atoms to be colored.
        :param excluded_index: the list of atom numbers of atoms to be excluded from coloring.
        :param bond_stereo: If true, add stereo information to bonds when constructing colors.
        :param resonance: If true, ignore the difference between double bonds and single bonds.
        :param backbone: If true, ignore bond types in the coloring.
        :return: None.
    """
        not_valid = self.invalid_symmetric_atoms(atoms_to_color, excluded_index, bond_stereo=bond_stereo,
                                                 resonance=resonance, backbone=backbone)
        while not_valid:
            atom_color_with_neighbors = self.generate_atom_color_with_neighbors(atoms_to_color, excluded=excluded_index,
                                                                                zero_core_color=False,
                                                                                zero_neighbor_color=False,
                                                                                resonance=resonance,
                                                                                bond_stereo=bond_stereo,
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

    def color_metal(self, bond_stereo: bool = False, resonance: bool = True, backbone: bool = False) -> None:
        """
        To color the metals in the compound. Here we just incorporate information of directly connected atoms.

        :param bond_stereo: If true, add bond stereo detail when constructing colors.
        :param resonance: If true, ignore difference between double and single bonds.
        :param backbone: If true, ignore the bond types.
        :return: None.
    """
        atom_color_with_neighbors = self.generate_atom_color_with_neighbors(self.metal_index, excluded=self.h_index,
                                                                            zero_core_color=True,
                                                                            zero_neighbor_color=False,
                                                                            resonance=resonance,
                                                                            bond_stereo=bond_stereo,
                                                                            backbone=backbone)
        for atom_index in self.metal_index:
            self.atoms[atom_index].color = atom_color_with_neighbors[atom_index]
      
    def color_h(self, bond_stereo: bool = False, resonance: bool = True, backbone: bool = False) -> None:
        """
        To color the H in the compound. Here we just incorporate information of directly connected atoms.

        :param bond_stereo:  If true, add bond stereo detail when constructing colors.
        :param resonance: If true, ignore difference between double and single bonds.
        :param backbone: If true, ignore bond types.
        :return: None.
    """
        atom_color_with_neighbors = self.generate_atom_color_with_neighbors(self.h_index, zero_core_color=True,
                                                                            zero_neighbor_color=False,
                                                                            resonance=resonance,
                                                                            bond_stereo=bond_stereo, backbone=backbone)
        for atom_index in self.h_index:
            self.atoms[atom_index].color = atom_color_with_neighbors[atom_index]

    def metal_color_identifier(self, details: bool = True) -> str:
        """
        To generate the metal coloring string representation.

        :param details: if true,  to use full metal color when constructing identifier.
        :return: the metal coloring string representation.
    """
        if details:  
            color_counter = collections.Counter([self.atoms[index].color for index in self.metal_index])
        else:
            color_counter = collections.Counter([self.atoms[index].color_0 for index in self.metal_index])
        return "".join(["({0})({1})".format(color_counter[key], key) for key in sorted(color_counter)])
    
    def h_color_identifier(self, details: bool = True) -> str:
        """
        To generate the H coloring string representation.

        :param details: if true, use the full H color when constructing identifier.
        :return: the H coloring string representation.
    """
        if details:
            color_counter = collections.Counter([self.atoms[index].color for index in self.h_index])
        else:
            color_counter = collections.Counter([self.atoms[index].color_0 for index in self.h_index])
        return "".join(["({0})({1})".format(color_counter[key], key) for key in sorted(color_counter)])

    def backbone_color_identifier(self, r_groups: bool = False) -> str:
        """
        To generate the backbone coloring string representation for this compound. Exclude Hs and metals.

        :param r_groups: whether to include the R group.
        :return: the coloring string representation for this compound.
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

    def get_chemical_details(self, excluded: list = None) -> list:
        """
        To get the chemical details of the compound, which include the atom stereo chemistry and bond stereo chemistry.
        This is to compare the compound with the same structures (or the same color identifiers).

        :param excluded: a list of atom indices to be ignored.
        :return: the list of chemical details in the compound.
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
    def compare_chemical_details(one_chemical_details: list, the_other_chemical_details: list) -> tuple:
        """
        To compare the chemical details of the two compounds.

        Then return the relationship between the two compounds.

        The relationship can be equivalent, generic-specific and loose, represented by 0, (-1, 1), 2

        :param one_chemical_details: the chemical details of one compound.
        :param the_other_chemical_details: the chemical details of the other compound.
        :return: the relationship between the two structures and the count of chemical details that cannot be mapped.
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
            # one_compound is more specific than the_other_compound
            return -1, len(one_more)
        elif the_other_more:
            return 1, len(the_other_more) # the_other_compound is more specific than one_compound

    def same_structure_relationship(self, the_other_compound) -> tuple:
        """
        To determine the relationship of two compounds with the same structure.

        :param the_other_compound: the other :class:`~md_harmonize.compound.Compound` entity.
        :return: the relationship and the atom mappings between the two compounds.
    """
        relationship, mis_count = self.compare_chemical_details(self.get_chemical_details(),
                                                                the_other_compound.get_chemical_details())
        # return relationship and atom mappings.
        return relationship, self.generate_atom_mapping_by_atom_color(the_other_compound)

    def generate_atom_mapping_by_atom_color(self, the_other_compound) -> dict:
        """
        To generate the atom mappings between the two compounds.

        Assume the two compounds have the same structure, so we can achieve atom mappings through atom colors.

        :param the_other_compound: the other :class:`~md_harmonize.compound.Compound` entity.
        :return: the atom mappings between the two compounds.
    """
        one_to_one_mappings = collections.defaultdict(list)
        for atom in self.heavy_atoms:
            for the_other_idx, the_other_atom in enumerate(the_other_compound.atoms):
                if atom.color == the_other_atom.color:
                    one_to_one_mappings[atom.atom_number].append(the_other_idx)
        return one_to_one_mappings

    def optimal_resonant_mapping(self, the_other_compound, mappings: list) -> tuple:
        """
        To find the optimal atom mappings for compound pairs that are resonant type.

        :param the_other_compound: the other :class:`~md_harmonize.compound.Compound` entity.
        :param mappings: the list of atom mappings between the two compounds detected by BASS.
        :return: the relationship and the atom mappings between the two compounds.
    """
        optimal_index = {1: -1, -1: -1, 2: -1, 0: -1}
        min_count = {1: float("inf"), -1: float("inf"), 2: float("inf"), 0: float("inf")}
        for i, mapping in enumerate(mappings):
            one_stereo_counts, the_other_stereo_counts = self.compare_chemical_details_with_mapping(the_other_compound,
                                                                                                    mapping)
            if not one_stereo_counts and not the_other_stereo_counts:
                optimal_index[0] = i
                min_count[0] = 0
            if not one_stereo_counts:
                if the_other_stereo_counts < min_count[1]:
                    optimal_index[-1] = i
                    min_count[1] = one_stereo_counts
            elif not the_other_stereo_counts:
                if one_stereo_counts < min_count[-1]:
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
    def determine_relationship(unmapped_count: dict) -> int:
        """
        To determine the relationship between two compounds when there are multiple possible atom mappings.

        We try to map as many details as possible.

        0: equivalent; 1: self is more generic than the other compound; -1: the other compound is more generic than self;
        2: either has chemical detail(s) that the other compound does not have.

        :param unmapped_count: the dictionary of relationship to the count of details that cannot be mapped.
        :return: the relationship between the two compounds.
    """
        if unmapped_count[0] == 0:
            # equivalent mapping.
            return 0
        if unmapped_count[-1] < unmapped_count[1]:
            return -1
        elif unmapped_count[-1] > unmapped_count[1]:
            return 1
        elif unmapped_count[-1] == unmapped_count[1] and unmapped_count[-1] != float("inf"):
            return 1
        return 2

    # def circular_pair_relationship(self, the_other_compound) -> tuple:
    #
    #     try:
    #         with tools.timeout(seconds=10):
    #             return self.circular_pair_relationship_helper(the_other_compound)
    #     except Exception as exception:
    #         return None, None

    def circular_pair_relationship(self, other_compound, seconds : int = 50) -> tuple:
        """
        To determine the relationship of two compounds with interchangeable circular and linear representations with time limit.

        :param other_compound: the other :class:`~md_harmonize.compound.Compound` entity.
        :param seconds: the timeout limit.
        :return: the relationship and the atom mappings between the two compounds.
        """

        try:
            relationship, mapping = tools.timeout(self.circular_pair_relationship_helper, (other_compound,), seconds=seconds)
            return relationship, mapping
        except Exception as exception:
            return None, None

    def circular_pair_relationship_helper(self, other_compound) -> tuple:
        """
        To determine the relationship of two compounds with interchangeable circular and linear representations.
        We first find the critical atoms that involve in the formation of the ring. There can be several possibilities.
        Then we break the ring, and restore the double bond in the aldehyde group that forms the ring.
        Finally, check if the updated structure is the same with the other compound. And determine the relationship
        between the two compounds as well as generate the atom mappings.

        :param other_compound: the other :class:`~md_harmonize.compound.Compound` entity.
        :return: the relationship and the atom mappings between the two compounds.
    """
        # default one compound should have a cycle.

        self.color_compound(r_groups=True, atom_stereo=False, bond_stereo=False)
        other_compound.color_compound(r_groups=True, atom_stereo=False, bond_stereo=False)

        optimal_mappings = {1: None, -1: None, 2: None, 0: None}
        min_count = {1: float("inf"), -1: float("inf"), 2: float("inf"), 0: float("inf")}
        # self.find_cycles()
        critical_atom_list = self.find_critical_atom_in_cycle()

        the_other_color = other_compound.backbone_color_identifier(r_groups=True) + \
                          other_compound.metal_color_identifier(details=False)
        for critical_atoms in critical_atom_list:
            self.break_cycle(critical_atoms)
            self.find_cycles()
            this_color = self.backbone_color_identifier(r_groups=True) + self.metal_color_identifier(details=False)
            if this_color == the_other_color:
                atom_mappings = self.generate_atom_mapping_by_atom_color(other_compound)
                excluded_atoms_the_other = list(set(itertools.chain.from_iterable([atom_mappings[i] for i in critical_atoms])))
                one_chemical_details = self.get_chemical_details(critical_atoms)
                the_other_chemical_details = other_compound.get_chemical_details(excluded_atoms_the_other)
                relationship, mis_count = self.compare_chemical_details(one_chemical_details, the_other_chemical_details)
                if mis_count < min_count[relationship]:
                    min_count[relationship] = mis_count
                    optimal_mappings[relationship] = atom_mappings
            self.restore_cycle(critical_atoms)
            self.find_cycles()
        if min(min_count.values()) < float("inf"):
            relationship = self.determine_relationship(min_count)
            return relationship, optimal_mappings[relationship]
        else:
            return None, None

    def break_cycle(self, critical_atoms: int) -> None:
        """
        To break the cycle caused by aldol reaction, which often occurs in the sugar.
        Two steps are involved: 1) remove the neighbors. 2) restore the double bond in the aldehyde group.

        :param critical_atoms: the three critical atoms that are involved in the ring formation.
        :return: None.
    """
        atom_o, atom_c, atom_oo = critical_atoms
        # break the cycle
        self.atoms[atom_o].remove_neighbors([atom_c])
        self.atoms[atom_c].remove_neighbors([atom_o])
        # update the single bond to double bond,
        self.bond_lookup[(atom_c, atom_oo)].update_bond_type("2")
        self.color_compound(r_groups=True, bond_stereo=False, atom_stereo=False, resonance=False,
                            isotope_resolved=False, charge=False)

    def restore_cycle(self, critical_atoms: list) -> None:
        """
        To restore the ring caused by aldol reaction.
        The reverse process of break_cycle.

        :param critical_atoms: the three atoms are involved in the aldol reaction.
        :return: None.
    """
        atom_o, atom_c, atom_oo = critical_atoms
        self.atoms[atom_o].add_neighbors([atom_c])
        self.atoms[atom_c].add_neighbors([atom_o])
        self.bond_lookup[(atom_c, atom_oo)].update_bond_type("1")
        self.color_compound(r_groups=True, bond_stereo=False, atom_stereo=False, resonance=False,
                            isotope_resolved=False,
                            charge=False)

    def find_critical_atom_in_cycle(self) -> list:
        """
        To find the C (atom_c) and O (atom_oo) in aldehyde group, as well as O (atom_o) in the hydroxy that are involved 
        in the ring formation. We need to break the bond between the atom_c and atom_o to form the linear transformation.
        Please check one example of aldol reaction in the sugar if the description is not confusing.

        :return: the list of critical atoms.
    """
        critical_atoms =[]
        for atom in self.atoms:
            # Two Os are connected to one atom and one O is in the cycle.
            if atom.in_cycle and atom.default_symbol == "O":
                for neighbor_index in atom.neighbors:
                    neighbor = self.atoms[neighbor_index]
                    for next_neighbor_index in neighbor.neighbors:
                        if next_neighbor_index != atom.atom_number and self.atoms[next_neighbor_index].default_symbol == "O" \
                                and self.bond_lookup[(next_neighbor_index, neighbor_index)].bond_type == "1":
                            critical_atoms.append([atom.atom_number, neighbor_index, next_neighbor_index])
        return critical_atoms

    def update_atom_symbol(self, index: list, updated_symbol: str) -> None:
        """
        To update the atom symbols. This is often used to remove/restore R group.

        :param index: the atom symbols of these indices to be updated.
        :param updated_symbol: the updated symbol.
        :return: None.
    """
        for i in index:
            self.atoms[i].update_symbol(updated_symbol)

    def validate_mapping_with_r(self, other_compound, one_rs: list, mapping: dict) -> bool:
        """
        To validate the atom mappings with r groups. 
        Here are two things we need to pay attention to:

        1) For the generic compound, the R group can be mapped to a branch or just H in the specific compound.

        2) For the specific compound, every unmatched branch needs to correspond to an R group in the generic compound.

        In other words, the generic compound can have extra R groups that have no matched branch, but the specific
        compound cannot have unmatched branches that don't correspond to any R groups.
        
        For the specific validation:

        1) We find all the linkages of R group and mapped atom in the compound, represented by the corresponding atom
        number in the other compound and the bond type. (We used the corresponding atom number in the other compound for
        the next comparison of the R linkages in the two compounds.

        2) For every mapped atom in the other compound, we need to find if it has neighbors that are not mapped. Then
        the atom should be linked to a R group. We represent the linkage by the atom number and the bond type.

        3) Based on the above validation criteria, we have to make sure that the R linkages in the other compound is the
        subset of the R linkages in this compound.
        
        :param other_compound: the other :class:`~md_harmonize.compound.Compound` entity.
        :param one_rs: the R groups in the compound.
        :param mapping: the atom mappings between the mapped parts of the two compounds.
        :return: bool whether the atom mappings are valid.
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
                one_r_linkages["{0}-{1}".format(mapping[neighbor_index], bond.bond_type)] += 1

        the_other_r_linkages = collections.Counter()
        for idx, atom in enumerate(other_compound.atoms):
            if idx in mapping.values():
                for neighbor_index in atom.neighbors:
                    if other_compound.atoms[neighbor_index].default_symbol != "H" and \
                            neighbor_index not in mapping.values():
                        bond = other_compound.bond_lookup[(idx, neighbor_index)]
                        the_other_r_linkages["{0}-{1}".format(idx, bond.bond_type)] += 1
        if all(one_r_linkages[x] >= the_other_r_linkages[x] for x in the_other_r_linkages):
            return True
        return False

    def compare_chemical_details_with_mapping(self, other_compound, mapping: dict) -> tuple:
        """
        To compare the chemical details of mapped atoms of the two compounds.
        This part targets compound pairs with resonance or r_group type.
        Only parts of chemicals need to be checked.
        1) atoms are not involved in resonance part or connected to R groups (both cases can be tested by the first
        layer atom coloring identifier).
        2) bond are formed by the atoms described above.

        :param other_compound: the other :class:`~md_harmonize.compound.Compound` entity.
        :param mapping: the mapped atoms between the two compounds.
        :return: the count of chemical details that cannot be mapped.
    """
        one_stereo_counts, the_other_stereo_counts = 0, 0
        one_consistent_atoms, the_other_consistent_atoms = set(), set()
        # check the atom stereochemistry
        for one_atom_index in mapping:
            one_atom = self.atoms[one_atom_index]
            the_other_atom = other_compound.atoms[mapping[one_atom_index]]
            if 1 in one_atom.color_layers and 1 in the_other_atom.color_layers and \
                    one_atom.color_layers[1] == the_other_atom.color_layers[1]:
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
                the_other_bond = other_compound.bond_lookup[(the_other_atom_1, the_other_atom_2)]
                if bond.bond_stereo == the_other_bond.bond_stereo:
                    continue
                if bond.bond_stereo != "0":
                    one_stereo_counts += 1
                if the_other_bond.bond_stereo != "0":
                    the_other_stereo_counts += 1
        return one_stereo_counts, the_other_stereo_counts

    def optimal_mapping_with_r(self, other_compound, one_rs: list, mappings: list) -> tuple:
        """
        To find the optimal mappings of compound pairs belonging to r_group type. In this case, multiple valid mappings
        can exist. We need to find the optimal one with the minimal unmapped chemical details.
        And the unmapped chemical details can exist in both compounds (generic or specific).
        The unmapped chemical details will determine the relationship of the compound pair.
        The priority: generic-specific, loose. 
        The relationship cannot be equivalent.

        :param other_compound: the other :class:`~md_harmonize.compound.Compound` entity.
        :param one_rs: the list of R groups in the compound.
        :param mappings: the atom mappings of the mapped parts in the two compounds.
        :return: the relationship and atom mappings between the two compounds.
    """
        optimal_index = {1: -1, -1: -1, 2: -1, 0: -1}
        min_count = {1: float("inf"), -1: float("inf"), 2: float("inf"), 0: float("inf")}
        for i, mm in enumerate(mappings):
            if self.validate_mapping_with_r(other_compound, one_rs, mm):
                one_stereo_counts, the_other_stereo_counts = self.compare_chemical_details_with_mapping(
                    other_compound, mm)
                if one_stereo_counts == 0 and the_other_stereo_counts < min_count[1]:
                    min_count[1] = the_other_stereo_counts
                    optimal_index[1] = i
                elif one_stereo_counts != 0 and one_stereo_counts + the_other_stereo_counts < min_count[2]:
                    min_count[2] = one_stereo_counts + the_other_stereo_counts
                    optimal_index[2] = i
        if min(min_count.values()) < float("inf"):
            relationship = self.determine_relationship(min_count)
            return relationship, mappings[optimal_index[relationship]]
        else:
            return None, None

    # def with_r_pair_relationship(self, the_other_compound) -> tuple:
    #
    #     try:
    #         with tools.timeout(seconds=10):
    #             return self.with_r_pair_relationship_helper(the_other_compound)
    #     except Exception as exception:
    #         return None, None

    def with_r_pair_relationship(self, other_compound, seconds: int = 50) -> tuple:
        """
        To find the relationship and the atom mappings between the two compounds that have r_groups type with a time limit.

        :param other_compound: the other :class:`~md_harmonize.compound.Compound` entity.
        :param seconds: the timeout limit.
        :return: the relationship and the atom mappings between the two compounds.
        """
        try:
            relationship, mapping = tools.timeout(self.with_r_pair_relationship_helper, (other_compound,), seconds=seconds)
            return relationship, mapping
        except Exception as exception:
            return None, None

    def with_r_pair_relationship_helper(self, other_compound) -> tuple:
        """
        To find the relationship and the atom mappings between the two compounds that have r_groups type.
        Several steps are involved:

        1) Ignore the R groups in the two compounds and find if one compound (generic compound) is included in the
        other compound (specific compound).

        2) If we can find the mappings, then we need to validate the mappings with the validate_mapping_with_r
        function.

        3）Then we get the optimal atom mappings of the mapped parts.

        4) We need to map the unmatched branches in the specific compound to the corresponding R group in the generic
        compound.

        :param other_compound: the other :class:`~md_harmonize.compound.Compound` entity.
        :return: the relationship and the atom mappings between the two compounds.
    """
        # self is the substructure, more generic, contain less chemical details. Apart from R groups, self is supposed
        # to contain smaller atoms, too.

        one_rs = list(self.r_groups)
        the_other_rs = list(other_compound.r_groups)
        self.update_atom_symbol(one_rs, "H")
        other_compound.update_atom_symbol(the_other_rs, "H")

        if len(self.heavy_atoms) > len(other_compound.heavy_atoms):
            return None, None

        if len(self.heavy_atoms) != 0:

            self.color_compound(r_groups=True, atom_stereo=False, bond_stereo=False)
            other_compound.color_compound(r_groups=True, atom_stereo=False, bond_stereo=False)
            one_to_one_mappings = self.find_mappings_reversed(other_compound, resonance=False, r_distance=True)

            # here we need to consider the r_distance atom color identifier, so we need to color compounds.
            relationship, optimal_mappings = self.optimal_mapping_with_r(other_compound, one_rs, one_to_one_mappings)
        else:
            relationship, optimal_mappings = 1, collections.defaultdict()
        if relationship:
            self.update_atom_symbol(one_rs, "R")
            other_compound.update_atom_symbol(the_other_rs, "R")
            return relationship, self.map_r_correspondents(one_rs, other_compound, optimal_mappings)

        return None, None

    def map_r_correspondents(self, one_rs: list, other_compound, mappings: dict) -> dict:
        # again, the self compound is substructure, we need to figure out the R group atom in self compound and its
        # corresponding atoms in the other compound.
        """ 
        To map the unmatched branches in the specific compound to the corresponding R group in the generic compound.

        :param one_rs: the list of R groups in the compound.
        :param other_compound: the other :class:`~md_harmonize.compound.Compound` entity.
        :param mappings: the atom mappings of the mapped parts in the two compounds.
        :return: the full atom mappings between the two compounds.
    """

        full_mappings = collections.defaultdict(list)

        if not mappings:
            r_index = one_rs[0]
            for atom in other_compound.heavy_atoms:
                full_mappings[r_index].append(atom.atom_number)
        else:
            for idx in one_rs:
                r_atom = self.atoms[idx]
                # it starts from R
                visited = set(mappings[neighbor_index] for neighbor_index in r_atom.neighbors)
                r_correspondents = [neighbor_index for index in visited for neighbor_index in
                                    other_compound.atoms[index].neighbors if neighbor_index not in mappings.values()]
                visited |= set(r_correspondents)
                while r_correspondents:
                    full_mappings[idx].extend(r_correspondents)
                    r_correspondents = other_compound.get_next_layer_neighbors(r_correspondents, visited)
            for idx in mappings:
                full_mappings[idx].append(mappings[idx])
        return full_mappings
