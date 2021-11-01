#!/usr/bin/python3

"""

KEGG_parser.py provides functions to parse KEGG data (including compound, reaction, kcf, and rclass).

"""

import collections
import glob
import re
import multiprocessing
import timeout_decorator

from . import compound
from . import reaction
from . import tools

def kegg_data_parser(data):
    """
    This is to parse KEGG data (reaction, rclass, compound) file to a dictionary.

    eg:
    ENTRY       R00259                      Reaction
    NAME        acetyl-CoA:L-glutamate N-acetyltransferase
    DEFINITION  Acetyl-CoA + L-Glutamate <=> CoA + N-Acetyl-L-glutamate
    EQUATION    C00024 + C00025 <=> C00010 + C00624
    RCLASS      RC00004  C00010_C00024
                RC00064  C00025_C00624
    ENZYME      2.3.1.1
    PATHWAY     rn00220  Arginine biosynthesis
                rn01100  Metabolic pathways
                rn01110  Biosynthesis of secondary metabolites
                rn01210  2-Oxocarboxylic acid metabolism
                rn01230  Biosynthesis of amino acids
    MODULE      M00028  Ornithine biosynthesis, glutamate => ornithine
                M00845  Arginine biosynthesis, glutamate => acetylcitrulline => arginine
    ORTHOLOGY   K00618  amino-acid N-acetyltransferase [EC:2.3.1.1]
                K00619  amino-acid N-acetyltransferase [EC:2.3.1.1]
                K00620  glutamate N-acetyltransferase / amino-acid N-acetyltransferase [EC:2.3.1.35 2.3.1.1]
                K11067  N-acetylglutamate synthase [EC:2.3.1.1]
                K14681  argininosuccinate lyase / amino-acid N-acetyltransferase [EC:4.3.2.1 2.3.1.1]
                K14682  amino-acid N-acetyltransferase [EC:2.3.1.1]
                K22476  N-acetylglutamate synthase [EC:2.3.1.1]
                K22477  N-acetylglutamate synthase [EC:2.3.1.1]
                K22478  bifunctional N-acetylglutamate synthase/kinase [EC:2.3.1.1 2.7.2.8]
    DBLINKS     RHEA: 24295
    ///

    :param data: the KEGG reaction description.
    :type data: :py:obj:`str`.
    :return: the dictionary of parsed KEGG data.
    :rtype: :py:obj:`dict`.
    """
    reaction_dict = collections.defaultdict(list)
    key = ""
    for line in data:
        if line.startswith("  "):
            reaction_dict[key].append(line[12:].rstrip())
        elif line.startswith("///"):
            break
        else:
            key = line[:12].strip()
            reaction_dict[key].append(line[12:].rstrip())
    return reaction_dict

def parse_equation(equation):
    """
    This is to parse the kegg reaction equation.

    eg:
    C00029 + C00001 + 2 C00003 <=> C00167 + 2 C00004 + 2 C00080

    :param equation: the equation string.
    :type equation: :py:obj:`str`.
    :return: the parsed KEGG reaction equation.
    :rtype: :py:obj:`dict`.
    """
    compound_pattern = "C....."
    one_side_coefficients = {}
    the_other_side_coefficients = {}
    current_side = one_side_coefficients
    last_compound = ""
    for token in equation.split():
        if re.search(compound_pattern, token):
            last_compound = token
            current_side[last_compound] = 1
        elif token.isdigit():
            current_side[last_compound] = int(token)
        elif token == "<=>":
            current_side = the_other_side_coefficients
    return one_side_coefficients, the_other_side_coefficients

def kegg_kcf_parser(kcf):
    """
    This is to parse KEGG kcf file to a dictionary.

    eg:
    ENTRY       C00013                      Compound
    ATOM        9
                1   P1b P    22.2269  -20.0662
                2   O2c O    23.5190  -20.0779
                3   O1c O    21.0165  -20.0779
                4   O1c O    22.2851  -21.4754
                5   O1c O    22.2617  -18.4642
                6   P1b P    24.8933  -20.0837
                7   O1c O    24.9401  -21.4811
                8   O1c O    26.1797  -20.0662
                9   O1c O    24.9107  -18.4582
    BOND        8
                1     1   2 1
                2     1   3 1
                3     1   4 1
                4     1   5 2
                5     2   6 1
                6     6   7 1
                7     6   8 1
                8     6   9 2
    ///

    :param kcf: the kcf text.
    :type kcf: :py:obj:`str`.
    :return: the dictionary of parsed kcf file.
    :rtype: :py:obj:`dict`.
    """
    compound_name, atom_count, atoms, bond_count, bonds = "", 0, [], 0, []

    for line in kcf:
        tokens = [item for item in line.split(' ') if item != '']
        if tokens[0] == "ENTRY":
            compound_name = tokens[1]
        elif tokens[0] == "ATOM":
            atom_count = int(tokens[1])
        elif tokens[0] == "BOND":
            bond_count = int(tokens[1])
        elif tokens[0] == "///":
            break
        else:
            if len(atoms) < atom_count:
                atoms.append({"atom_symbol": tokens[2], "atom_number": int(tokens[0])-1, "x": tokens[3], "y": tokens[4], "kat": tokens[1]})
            elif len(bonds) < bond_count:
                bonds.append({"first_atom_number":tokens[1], "second_atom_number": tokens[2], "bond_type": tokens[3]})
    return {"compound_name": compound_name, "atoms": atoms, "bonds": bonds}


# Define the reaction center in the rclass definition.
reaction_center = collections.namedtuple('reaction_center', ['i', 'kat', 'label', 'match', 'difference'])

class RpairParser:

    """This is to get one to one atom mappings between two compounds based on the rclass definition.

        Several steps are involved in this process:
        1) The rclass definition can have several pieces. Each piece describes a center atom (R) and its connected atoms.
        The connected atoms can stay the same (M) or change (D) between the two compound structures.
        2) First we need to find the center atoms based on the rlcass descriptions.
        3) For each center atom, there are can multiple candidates. In other words, based on the RDM description,
        a bunch of atoms in the compound can meet the descriptions. (One simple case are the symmetric compounds).
        4) Therefore, we need to generate the all the combinations for the center atoms in a compound.
        eg: if there are three atom centers, each center has several candidates:
            center 1: [0, 1, 2]; center 2: [5, 6]; center 3: [10, 11]
            The combinations for the center atoms:
            [0, 5, 10], [0, 5, 11], [0, 6, 10], [0, 6, 11], [1, 5, 10], [1, 5, 11], [1, 6, 10], [1, 6, 11], [2, 5, 10],
            [2, 5, 11], [2, 6, 10], [2, 6, 11]
        5) Next, we need to find the one to one atom mappings between the two compounds based on the mapped center atoms.
        6) To solve this issue, we first disassemble each compound into several different components. This is due to the
        difference atoms in the two compounds, i.e. broken bonds.
        7) Then we need to find the mappings between each disassembled component, and concatenate the mappings of all the components.
        8) To find the one to one atom mappings, we use the BASS algorithm. We assume the mapped component have the same
        structure since we have already removed the different parts. However, here we only map the backbone of the structure (in other words,
        we simply all the bond type to 1) due to bond change (double bond to single bond or triple bond to single bond)
        9) To ensure the optimal mappings, we count the mapped atoms with changed local environment and choose the mapping
        with minimal changed local colors.
    """

    def __init__(self, rclass_name, rclass_definitions, one_compound, the_other_compound):
        """
        RpairParser initializer.

        :param rclass_name: the rclass name.
        :type rclass_name: :py:obj:`str`.
        :param rclass_definitions: a list of rclass definitions.
        :type rclass_definitions: :py:obj:`list`.
        :param one_compound: one compound involved in the pair.
        :type one_compound: :class:`~MDH.compound.Compound`.
        :param the_other_compound: the other compound involved in the pair.
        :type the_other_compound: :class:`~MDH.compound.Compound`.
        """
        self.rclass_name = rclass_name
        self.rclass_definitions = rclass_definitions
        self.one_compound = one_compound
        self.the_other_compound = the_other_compound

    @staticmethod
    def generate_kat_neighbors(compound):
        """
        To generate the atom neighbors represented by KEGG atom type for each atom in the compound.
        This is used to find the center atom.
        We used KEGG atom type since the descriptions of atoms in rclass definitions use KEGG atom type.

        :param compound: the compound entity.
        :type compound: :class:`~MDH.compound.Compound`.
        :return: the list of atom with its neighbors.
        :rtype: :py:obj:`list`.
        """
        atoms = []
        for atom in compound.atoms:
            kat_neighbors = [compound.atoms[i].kat for i in atom.neighbors if compound.atoms[i].default_symbol != "H"]
            kat_neighbors.sort()
            atoms.append((atom.kat, kat_neighbors))
        return atoms

    @staticmethod
    def search_target_atom(atoms, target):
        """
        To search the target atom from a list of atoms.

        :param atoms: a list of atoms to be searched.
        :type atoms: :py:obj:`list`.
        :param target: the target atom to be searched.
        :type target: :py:obj:`set`.
        :return: the list of atom number that match the target.
        :rtype: :py:obj:`list`.
        """
        target_index = []
        for i, atom in enumerate(atoms):
            if atom == target:
                target_index.append(i)
        return target_index

    @staticmethod
    def create_reaction_centers(i, kat, difference, the_other_difference, match, the_other_match):
        """
        To create the center atom based on its connected atoms and the its counterpart atom in the other compound.

        :param i: the ith rclass definition.
        :type i: :py:obj:`int`.
        :param kat: KEGG atom type of the center atom.
        :type kat: :py:obj:`str`.
        :param difference: the list of KEGG atom type of different connected atoms
        :type difference: :py:obj:`list`.
        :param the_other_difference: the list of KEGG atom type of different connected atoms of the other compound.
        :type the_other_difference: :py:obj:`list`.
        :param match: the list of KEGG atom type of the matched connected atoms.
        :type match: :py:obj:`list`.
        :param the_other_match: the list of KEGG atom type of the matched connected atoms of the other compound.
        :type the_other_match: :py:obj:`list`.
        :return: the constructed reaction center.
        :rtype: py:class:`~collections.namedtuple`.
        """
        match = list(zip(match, the_other_match))
        difference = list(zip(difference, the_other_difference))
        neighbors = [item for item in match if item != "*"]
        neighbors.sort()
        label = (kat, tuple(neighbors))
        return reaction_center(i, kat, label, match, difference)

    def find_center_atoms(self):
        """
        Example of rclass definition:
        C8x-C8y:*-C1c:N5y+S2x-N5y+S2x
        The RDM pattern is defined as KEGG atom type changes at the reaction center (R), the difference region (D),
        and the matched region (M) for each reactant pair. It characterizes chemical structure transformation patterns
        associated with enzymatic reactions.

        :return: the list of reaction centers and their corresponding candidate atoms.
        :rtype: :py:obj:`list`.
        """
        left_reaction_centers, right_reaction_centers = [], []
        left_center_candidates, right_center_candidates = [], []
        left_kat_neighbors = self.generate_kat_neighbors(self.one_compound)
        right_kat_neighbors = self.generate_kat_neighbors(self.the_other_compound)
        for i, definition in enumerate(self.rclass_definitions):
            center, difference, match = definition.split(":")
            left_center, right_center = center.split("-")
            left_difference = [ item for item in difference.split("-")[0].split("+") ]
            right_difference = [ item for item in difference.split("-")[1].split("+") ]
            left_match = [ item for item in match.split("-")[0].split("+") ]
            right_match = [ item for item in match.split("-")[1].split("+") ]
            left_neighbors = [item for item in left_difference if item != "*"] + [item for item in left_match if item != "*"]
            right_neighbors = [item for item in right_difference if item != "*"] + [item for item in right_match if item != "*"]
            left_neighbors.sort()
            right_neighbors.sort()
            left_atom_index = self.search_target_atom(left_kat_neighbors, (left_center, left_neighbors))
            right_atom_index = self.search_target_atom(right_kat_neighbors, (right_center, right_neighbors))
            left_reaction_centers.append(self.create_reaction_centers(i, left_center, left_difference, right_difference, left_match, right_match))
            right_reaction_centers.append(self.create_reaction_centers(i, right_center, right_difference, left_difference, right_match, left_match))
            left_center_candidates.append(left_atom_index)
            right_center_candidates.append(right_atom_index)
        return left_center_candidates, right_center_candidates, left_reaction_centers, right_reaction_centers

    @staticmethod
    def get_center_list(center_atom_index):
        """
        To generate all the combinations of reaction centers.
        
        :param center_atom_index: list of atom index list for each reaction centers. eg: three reaction centers:
        [[0, 1, 2], [5, 6], [10, 11]].
        :type center_atom_index: :py:obj:`list`.
        :return: the list of combined reaction centers. eg: [[0, 5, 10], [0, 5, 11], [0, 6, 10], [0, 6, 11], [1, 5, 10],
        [1, 5, 11], [1, 6, 10], [1, 6, 11], [2, 5, 10], [2, 5, 11], [2, 6, 10], [2, 6, 11]]
        :rtype: :py:obj:`list`.
        """
        combines = []
        def dfs(i, seen):
            if i == len(center_atom_index):
                combines.append(list(seen))
                return
            for idx in center_atom_index[i]:
                if idx not in seen:
                    seen.append(idx)
                    dfs(i+1, seen)
                    seen.pop()
        dfs(0, [])
        return combines

    @staticmethod
    def remove_different_bonds(compound, center_atom_numbers, reaction_centers):
        """
        To remove the bonds connecting to different atoms. For each reaction center, multiple atoms can be the different atoms.
        We need to get all the combinations.
        
        :param compound: the :class:`~MDH.compound.Compound` entity.
        :type compound: :class:`~MDH.compound.Compound`.
        :param center_atom_numbers: the list of atom numbers of center atom in the compound.
        :type center_atom_numbers: :py:obj:`list`.
        :param reaction_centers: the list of reaction center descriptions of the compound.
        :type reaction_centers: :py:obj:`list`.
        :return: the list of bonds (represented by the atom numbers that forming the bond) that needs to be removed
        based on the RDM descriptions.
        :rtype: :py:obj:`list`.
        """
        removed_bonds = [[]]
        for i, idx in enumerate(center_atom_numbers):
            for different_atom in reaction_centers[i].difference:
                if different_atom[0] != "*":
                    removed_bonds_update = []
                    possible_bonds = []
                    for neighbor_index in compound.atoms[idx].neighbors:
                        if compound.atoms[neighbor_index].kat == different_atom[0]:
                            possible_bonds.append((idx, neighbor_index))
                        for possible_bond in possible_bonds:
                            for pre_removed in removed_bonds:
                                removed_bonds_update.append(pre_removed + [possible_bond])
                    removed_bonds = removed_bonds_update
        return removed_bonds

    def generate_atom_mappings(self):
        """
        To generate the one to one atom mappings of the compound pair.

        :return: the list of atom mappings.
        :rtype: :py:obj:`list`.
        """
        left_center_candidates, right_center_candidates, left_reaction_centers, right_reaction_centers = self.find_center_atoms()
        left_centers_list = self.get_center_list(left_center_candidates)
        right_center_list = self.get_center_list(right_center_candidates)

        minimum_miss_count = float("inf")
        optimal_atom_mappings = {}
        for left_centers in left_centers_list:
            for right_centers in right_center_list:
                # After removing the bonds, the compound will be disassembled into several parts. Try to match these pieces.
                left_removed_bond_list = self.remove_different_bonds(self.one_compound, left_centers, left_reaction_centers)
                right_removed_bonds_list = self.remove_different_bonds(self.the_other_compound, right_centers, right_reaction_centers)
                for left_removed_bonds in left_removed_bond_list:
                    for right_removed_bonds in right_removed_bonds_list:
                        this_mapping = self.map_components(left_removed_bonds, right_removed_bonds,
                                                           left_centers, right_centers)
                        if this_mapping:
                            miss_count = self.count_different_atom_identifiers(this_mapping)
                            if len(optimal_atom_mappings) < len(this_mapping) or \
                                    (len(optimal_atom_mappings) == len(this_mapping) and minimum_miss_count > miss_count):
                                minimum_miss_count = miss_count
                                optimal_atom_mappings = this_mapping

        one_to_one_mappings = []
        for from_idx in optimal_atom_mappings:
            to_idx = optimal_atom_mappings[from_idx]
            one_to_one_mappings.append(((self.one_compound.compound_name, from_idx), (self.the_other_compound.compound_name, to_idx)))
        return one_to_one_mappings

    @staticmethod
    def detect_components(compound, removed_bonds, center_atom_index):
        """

        :param compound:
        :param removed_bonds:
        :param center_atom_index:
        :return:
        """
        graph = collections.defaultdict(list)
        for bond in compound.bonds:
            i, j = bond.first_atom_number, bond.second_atom_number
            if (i, j) not in removed_bonds and (j, i) not in removed_bonds:
                if compound.atoms[i].default_symbol != "H" and compound.atoms[j].default_symbol != "H":
                    graph[i].append(j)
                    graph[j].append(i)

        nodes = set(list(graph.keys()) + center_atom_index)
        visited = set()
        components = []
        for node in nodes:
            if node not in visited:
                ques = collections.deque([node])
                component = []
                visited.add(node)
                while ques:
                    node = ques.popleft()
                    component.append(node)
                    for neighbor_node in graph[node]:
                        if neighbor_node not in visited:
                            visited.add(neighbor_node)
                            ques.append(neighbor_node)
                components.append(component)
        return components

    @staticmethod
    def get_component_pairs(left_components, right_components):
        """
        To get the component pairs with the same number of atoms.
        :param left_components:
        :param right_components:
        :return:
        """
        component_pairs = []
        for left_component in left_components:
            for right_component in right_components:
                if len(left_component) == len(right_component):
                    component_pairs.append((left_component, right_component))
        return component_pairs

    @staticmethod
    def construct_component(cpd, atom_index, removed_bonds):
        """
        To construct a partial compound based on the atom index and the removed bonds.

        :param cpd:
        :type cpd:
        :param atom_index:
        :param removed_bonds:
        :return:
        """
        atoms = [cpd.atoms[idx].clone() for idx in atom_index]
        idx_dict = {int(atom.atom_number): i for i, atom in enumerate(atoms)}
        for i, atom in enumerate(atoms):
            atom.update_atom_number(i)
        bonds = []
        for bond in cpd.bonds:
            atom_1, atom_2 = bond.first_atom_number, bond.second_atom_number
            if atom_1 in atom_index and atom_2 in atom_index and (atom_1, atom_2) not in removed_bonds and \
                    (atom_2, atom_1) not in removed_bonds:
                cloned_bond = bond.clone()
                cloned_bond.update_first_atom(idx_dict[atom_1])
                cloned_bond.update_second_atom(idx_dict[atom_2])
                bonds.append(cloned_bond)
        return compound.Compound("partial_compound", atoms, bonds)

    @staticmethod
    def preliminary_atom_mappings_check(left_partial_compound, right_partial_compound):
        """
        To roughly evaluate if the atoms between the two partial compounds can be mapped.
        :param left_partial_compound:
        :param right_partial_compound:
        :return:
        """
        left_partial_compound.color_compound(r_groups=False, bond_stereo=False, atom_stereo=False, resonance=True)
        right_partial_compound.color_compound(r_groups=False, bond_stereo=False, atom_stereo=False, resonance=True)
        left, mapped = 0, 0
        for atom_left in left_partial_compound.atoms:
            if atom_left.default_symbol != "H":
                left += 1
                for atom_right in right_partial_compound.atoms:
                    if atom_left.color_layers == atom_right.color_layers:
                        mapped += 1
                        break
        #print(left, mapped)
        return left == mapped

    def map_components(self, left_removed_bonds, right_removed_bonds, left_index_comb, right_index_comb):
        """

        :param left_removed_bonds:
        :param right_removed_bonds:
        :param left_index_comb:
        :param right_index_comb:
        :return:
        """

        left_components = self.detect_components(self.one_compound, left_removed_bonds, left_index_comb)
        right_components = self.detect_components(self.the_other_compound, right_removed_bonds, right_index_comb)
        component_pairs = self.get_component_pairs(left_components, right_components)
        atom_mappings = []
        for left_component, right_component in component_pairs:
            left_partial_compound = self.construct_component(self.one_compound, left_component, left_removed_bonds)
            right_partial_compound = self.construct_component(self.the_other_compound, right_component, right_removed_bonds)
            if not self.preliminary_atom_mappings_check(left_partial_compound, right_partial_compound):
                continue
            mappings = left_partial_compound.find_mappings(right_partial_compound, resonance=True, r_distance=False)
            one_to_one_mappings_list = left_partial_compound.generate_one_to_one_mappings(right_partial_compound, mappings)
            optimal_one_to_one_mappings = None
            minimum_miss_count = float("inf")
            for one_to_one_mappings in one_to_one_mappings_list:
                original_index_mappings = {}
                for i in one_to_one_mappings:
                    original_index_mappings[left_component[i]] = right_component[one_to_one_mappings[i]]
                    if self.validate_component_atom_mappings(left_index_comb, right_index_comb, original_index_mappings):
                        miss_count = self.count_different_atom_identifiers(original_index_mappings)
                        if miss_count < minimum_miss_count:
                            minimum_miss_count = miss_count
                            optimal_one_to_one_mappings = original_index_mappings
            if optimal_one_to_one_mappings:
                atom_mappings.append(optimal_one_to_one_mappings)
        return self.combine_component_atom_mappings(atom_mappings)

    def combine_component_atom_mappings(self, atom_mappings):
        """

        :param atom_mappings:
        :return:
        """
        groups = collections.defaultdict(list)
        collective_mappings = {}
        for i, atom_mapping in enumerate(atom_mappings):
            groups[len(atom_mapping)].append(i)
        for g in groups:
            if len(groups[g]) == 1:
                collective_mappings.update(atom_mappings[groups[g][0]])
            else:
                orders = groups[g]
                orders.sort(key=lambda x: self.count_different_atom_identifiers(atom_mappings[x]))
                for idx in orders:
                    if all(atom_index not in collective_mappings for atom_index in atom_mappings[idx]) and \
                            all(atom_index not in collective_mappings.values() for atom_index in atom_mappings[idx].values()):
                        collective_mappings.update(atom_mappings[idx])
        return collective_mappings

    @staticmethod
    def validate_component_atom_mappings(left_centers, right_centers, component_atom_mappings):
        """
        To check if the mapped the atoms can corresponds to the mapped reaction center atoms.
        
        :param left_centers: the list of center atom index in the left compound.
        :type left_centers: :py:obj:`list`.
        :param right_centers: the list of center atom index in the right compound.
        :type right_centers: :py:obj:`list`.
        :param component_atom_mappings: the one to one atom mappings of one component.
        :type component_atom_mappings: :py:obj:`dict`.
        :return: bool whether the mappings are valid.
        :rtype: :py:class:`bool`.
        """
        for i, j in zip(left_centers, right_centers):
            if i in component_atom_mappings and j != component_atom_mappings[i]:
                return False
        return True

    def count_different_atom_identifiers(self, one_to_one_mappings):
        """
        To count the mapped atoms with changed local atom identifier. The different atoms (D in RCLASS definitions)
        can cause change of local environment, which can change the atom identifier.

        :param one_to_one_mappings: the dictionary of atom mappings between the two compounds.
        :type one_to_one_mappings: :py:obj:`dict`.
        :return: the total number of mapped atoms with different local identifier.
        :rtype: :py:obj:`int`.
        """
        count = 0
        for idx_1 in one_to_one_mappings:
            idx_2 = one_to_one_mappings[idx_1]
            if self.one_compound.atoms[idx_1].color_layers and self.the_other_compound.atoms[idx_2].color_layers and \
                    self.one_compound.atoms[idx_1].color_layers[1] == self.the_other_compound.atoms[idx_2].color_layers[1]:
                continue
            count += 1
        return count


def create_compound_kcf(kcf_file):
    """
    To construct compound entity based on the KEGG kcf file.

    :param kcf_file: the filename contains kcf text.
    :type kcf_file: :py:obj:`str`.
    :return: the constructed compound entity.
    :rtype: :class:`~MDH.compound.Compound`.
    """
    kcf_dict = kegg_kcf_parser(tools.open_text(kcf_file).split("\n"))
    atoms = [compound.Atom(atom["atom_symbol"], atom["atom_number"], x=atom["x"], y=atom["y"], kat=atom["kat"]) for
             atom in kcf_dict["atoms"]]
    bonds = [compound.Bond(bond["first_atom_number"], bond["second_atom_number"], bond["bond_type"]) for bond in kcf_dict["bonds"]]
    
    return compound.Compound(kcf_dict["compound_name"], atoms, bonds)

def add_kat(compound, kcf_compound):
    """
    To add kegg atom type to each atom in the compound based on the corresponding kcf compound.

    :param compound: a compound entity.
    :type compound: :class:`~MDH.compound.Compound`.
    :param kcf_compound: the corresponding compound entity parsed from KCF file.
    :type kcf_compound: :class:`~MDH.compound.Compound`.
    :return: None.
    :rtype: :py:obj:`None`.
    """
    compound.color_compound(r_groups=True, bond_stereo=False, atom_stereo=False, resonance=False)
    kcf_compound.color_compound(r_groups=True, bond_stereo=False, atom_stereo=False, resonance=False)
    color_kat = { atom.color: atom.kat for atom in kcf_compound.atoms if atom.default_symbol != "H" }

    for atom in compound.atoms:
        if atom.default_symbol != "H":
            if atom.color not in color_kat:
                print(atom.atom_number)
            atom.kat = color_kat[atom.color]

# when we create the kegg reaction, we need to parse the atom mappings based on rclass!
# To avoid parsing the same rclass repeatedly, let's parse the rclass first, and look it up when we need.
def create_reactions(reaction_directory, compounds, atom_mappings):
    """
    To create KEGG `~MDH.reaction.Reaction` entities.

    :param reaction_directory: the directory that stores all the reaction files.
    :type reaction_directory: :py:obj:`str`.
    :param compounds: a dictionary of :class:`~MDH.compound.Compound` entities.
    :type compounds: :py:obj:`dict`.
    :param atom_mappings: the compound pair name and the its atom mappings.
    :type atom_mappings: :py:obj:`dict`.
    :return: the constructed `~MDH.reaction.Reaction` entities.
    :rtype: :py:obj:`list`.
    """
    # here we compounds, rlcass descriptions, and reactions.
    reaction_files = glob.glob(reaction_directory+"*")
    reactions = []
    for reaction_file in reaction_files:
        this_reaction = kegg_data_parser(tools.open_text(reaction_file).split("\n"))
        reaction_name = this_reaction["ENTRY"][0].split()[0]

        raw_ecs = this_reaction["ENZYME"] if "ENZYME" in this_reaction else []
        ecs = collections.defaultdict(list)
        for line in raw_ecs:
            ec_numbers = line.split()
            for ec in ec_numbers:
                numbers = ec.split(".")
                # some numbers in the ec are not specified. like "3.5.99.-"
                while numbers and numbers[-1] == "-":
                    numbers.pop()
                ecs[len(numbers)].append(".".join(numbers))
        if not ecs:
            # ORTHOLOGY can contain information of ecs. Please check the example of reaction.
            for line in this_reaction["ORTHOLOGY"]:
                ec_pattern = "\[EC:.*\]"
                ec_numbers = re.findall(ec_pattern, line)[0][4:-1].split()
                for ec in ec_numbers:
                    numbers = ec.split(".")
                    while numbers and numbers[-1] == "-":
                        numbers.pop()
                    ecs[len(numbers)].append(".".join(numbers))

        one_side_coefficients, the_other_side_coefficients = parse_equation(this_reaction["EQUATION"][0])
        one_side_compounds = [compounds["cpd:" + compound_name] for compound_name in one_side_coefficients]
        the_other_side_compounds = [compounds["cpd:" + compound_name] for compound_name in the_other_side_coefficients]
        one_side_coefficients.update(the_other_side_coefficients)

        this_atom_mappings = []
        for line in this_reaction["RCLASS"]:
            rclass, rpair = line.split()
            cpd_1, cpd_2 = rpair.split("_")
            key_1 = "_".join([rclass, cpd_1, cpd_2])
            key_2 = "_".join([rclass, cpd_2, cpd_1])
            if key_1 in atom_mappings:
                this_atom_mappings.extend(atom_mappings[key_1])
            elif key_2 in atom_mappings:
                this_atom_mappings.extend(atom_mappings[key_2])
        reactions.append(reaction.Reaction(reaction_name, one_side_compounds, the_other_side_compounds, ecs,
                                           this_atom_mappings, one_side_coefficients))
    return reactions

@timeout_decorator.timeout(200)
def compound_pair_mappings(rclass_name, rclass_definitions, one_compound, the_other_compound):
    """
    To get the atom mappings between two compounds based on the rclass definitions.

    :param rclass_name: the name of the rclass.
    :type rclass_name: :py:obj:`str`.
    :param rclass_definitions: the list of rclass definitions.
    :type rclass_definitions: :py:obj:`list`.
    :param one_compound: one compound entity involved in the compound pair.
    :type one_compound: :class:`~MDH.compound.Compound`.
    :param the_other_compound: the other compound entity involved in the compound pair.
    :type the_other_compound: :class:`~MDH.compound.Compound`.
    :return: the compound pair name and the its atom mappings.
    :rtype: :py:obj:`str` and :py:obj:`list`.
    """
    atom_mappings = []
    try:
        one_mappings = RpairParser(rclass_name, rclass_definitions, one_compound, the_other_compound).generate_atom_mappings()
        the_other_mappings = RpairParser(rclass_name, rclass_definitions, the_other_compound, one_compound).generate_atom_mappings()
        atom_mappings = one_mappings if len(one_mappings) > len(the_other_mappings) else the_other_mappings
    except Exception:
        print(rclass_name + "_" + one_compound.compound_name[4:] + "_" + the_other_compound.compound_name[4:] + "can hardly be parsed.")
        pass
    return one_compound.compound_name[4:] + "_" + the_other_compound.compound_name[4:], atom_mappings

def create_atom_mappings(rclass_directory, compounds):
    """
    To generate the atom mappings between compounds based on RCLASS definitions.

    :param rclass_directory: the directory that stores the rclass files.
    :type rclass_directory: :py:obj:`str`.
    :param compounds: a dictionary of :class:`~MDH.compound.Compound` entities.
    :type compounds: :py:obj:`dict`.
    :return: the atom mappings of compound pairs.
    :rtype: :py:obj:`dict`.
    """
    rclass_files = glob.glob(rclass_directory + "*")
    atom_mappings = collections.defaultdict(dict)
    for rclass_file in rclass_files:
        this_rclass = kegg_data_parser(tools.open_text(rclass_file).split("\n"))
        rclass_definitions = this_rclass["DEFINITION"]
        rclass_name = this_rclass["ENTRY"][0].split()[0]
        compound_pairs = []
        for line in this_rclass["RPAIR"]:
            tokens = line.split()
            for token in tokens:
                one_compound_name, the_other_compound_name = token.split("_")
                if "cpd:" + one_compound_name in compounds and "cpd:" + the_other_compound_name in compounds:
                    compound_pairs.append((compounds["cpd:" + one_compound_name], compounds["cpd:" + the_other_compound_name]))

        with multiprocessing.Pool() as pool:
            results = pool.starmap(compound_pair_mappings, ((rclass_name, rclass_definitions, one_compound, the_other_compound)
                                                            for one_compound, the_other_compound in compound_pairs))
        for name, mapping in results:
            atom_mappings[rclass_name + "_" + name] = mapping
            if not mapping:
                print("empty mappings", rclass_name + "_" + name)
    return atom_mappings

