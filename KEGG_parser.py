#!/usr/bin/python3

import collections
import compound
import glob
import tools
from pathlib import Path
import ctfile
import openbabel
import aromatize
import re
import reaction

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

    :param str data: the KEGG reaction description.
    :return: the dictionary of the parsed KEGG reaction.
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
    :return:
    """

    one_side, the_other_side = equation.split("<=>")
    compound_pattern = "C....."
    one_side_coefficients = {}
    last_compound = ""
    for token in one_side.split():
        if re.search(compound_pattern, token):
            last_compound = token
            one_side_coefficients[last_compound] = 1
        elif token.isdigit():
            one_side_coefficients[last_compound] = int(token)
    the_other_side_coefficients = {}
    last_compound = ""
    for token in the_other_side.split():
        if re.search(compound_pattern, token):
            last_compound = token
            the_other_side_coefficients[last_compound] = 1
        elif token.isdigit():
            the_other_side_coefficients[last_compound] = int(token)

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

    :param kcf: the kcf text
    :return: the dictionary of parsed kcf file.
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
            else:
                bonds.append({"first_atom_number":tokens[1], "second_atom_number": tokens[2], "bond_type": tokens[3]})
    return {"compound_name": compound_name, "atoms": atoms, "bonds": bonds}


reaction_center = collections.namedtuple('reaction_center', ['i', 'kat', 'label', 'match', 'difference'])

class RpairParser:

    """This is to get one to one atom mappings between two compounds based on the rclass definition.

        Several steps are involved in this process:
        1) First we need to find the center atoms based on the rlcass descriptions. The number of center atoms equals to
        the number of rclass RMD descriptions.
        2) For each center atom, there are can multiple candidate atoms. In other words, based on the RDM description,
        several different atoms in the compound can meet the descriptions.
        3) Therefore, we need to generate the all the combinations for the center atoms in a compound. In this case, each
        compound has a set of list of center atom index.
        4) Next, we need to find the one to one atom mappings between the two compounds based on the mapped center atoms.
        5) To solve this issue, we first disassemble each compound into several different components. This is due to the
        difference atoms in the two compounds.

        Here, we can see that several different combinations are involved. Please pay attention to this!
    """

    def __init__(self, rclass_name, rclass_definitions, one_compound, the_other_compound):
        """

        :param rclass_name:
        :param rclass_definitions:
        :param one_compound:
        :param the_other_compound:
        """
        self.rclass_name = rclass_name
        self.rclass_definitions = rclass_definitions
        self.one_compound = one_compound
        self.the_other_compound = the_other_compound

    @staticmethod
    def generate_kat_neighbors(compound):
        """
        To generate the
        :param compound:
        :return:
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

        :param atoms:
        :param target:
        :return:
        """
        target_index = []
        for i, atom in enumerate(atoms):
            if atom == target:
                target_index.append(i)
        return target_index

    @staticmethod
    def create_reaction_center(i, kat, difference, the_other_difference, match, the_other_match):
        """

        :param i:
        :param kat:
        :param difference:
        :param the_other_difference:
        :param match:
        :param the_other_match:
        :return:
        """
        match = list(zip(match, the_other_match))
        difference = list(zip(difference, the_other_difference))
        neighbors = [item for item in match if item != "*"]
        neighbors.sort()
        label = (kat, tuple(neighbors))
        return reaction_center(i, kat, label, match, difference)

    @staticmethod
    def center_index_combinations(center_atom_index):
        """
        To generate the combinations of reaction center index.
        :param center_atom_index: list of atom index list for each reaction centers.
        :return:
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
    def remove_difference_bonds(compound, reaction_center_index, reaction_centers):
        """
        To remove the bonds connecting to different atoms. For each reaction center, multiple atoms can be the different atoms.
        We need to get all the combinations.
        :param compound:
        :param reaction_center_index:
        :param reaction_centers:
        :return:
        """
        removed_bonds = [[]]
        for i, idx in enumerate(reaction_center_index):
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

    def find_center_atoms(self):
        """

        Example of rclass definition:
        C8x-C8y:*-C1c:N5y+S2x-N5y+S2x
        The RDM pattern is defined as KEGG atom type changes at the reaction center (R), the difference region (D),
        and the matched region (M) for each reactant pair. It characterizes chemical structure transformation patterns
        associated with enzymatic reactions.

        :return:
        """
        left_center_atoms, right_center_atoms = [], []
        left_center_atom_index, right_center_atom_index = [], []
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
            left_center_atoms.append(self.create_reaction_center(i, left_center, left_difference, right_difference, left_match, right_match))
            right_center_atoms.append(self.create_reaction_center(i, right_center, right_difference, left_difference, right_match, left_match))
            left_center_atom_index.append(left_atom_index)
            right_center_atom_index.append(right_atom_index)
        return left_center_atom_index, right_center_atom_index, left_center_atoms, right_center_atoms

    def map_center_atoms(self):
        """

        :return:
        """

        left_center_atom_index, right_center_atom_index, left_center_atoms, right_center_atoms = self.find_center_atoms()
        left_index_combs = self.center_index_combinations(left_center_atom_index)
        right_index_combs = self.center_index_combinations(right_center_atom_index)

        minimum_miss_count = float("inf")
        optimal_atom_mappings = {}
        for left_index_comb in left_index_combs:
            for right_index_comb in right_index_combs:
                # After removing the bonds, the compound will be disassembled into several parts. Try to match these pieces.
                left_removed_bond_list = self.remove_difference_bonds(self.one_compound, left_index_comb, left_center_atoms)
                right_removed_bonds_list = self.remove_difference_bonds(self.the_other_compound, right_index_comb, right_center_atoms)
                for left_removed_bonds in left_removed_bond_list:
                    for right_removed_bonds in right_removed_bonds_list:
                        this_mapping = self.map_separate_components(left_removed_bonds, right_removed_bonds,
                                                                    left_index_comb, right_index_comb)
                        if this_mapping:
                            miss_count = self.count_different_atom_identifiers(this_mapping)
                            if len(optimal_atom_mappings) < len(this_mapping) or \
                                    (len(optimal_atom_mappings) == len(this_mapping) and minimum_miss_count > miss_count):
                                minimum_miss_count = miss_count
                                optimal_atom_mappings = this_mapping
        return self.generate_one_to_one_mappings(optimal_atom_mappings)

    def generate_one_to_one_mappings(self, atom_mappings):
        """

        :param atom_mappings:
        :return:
        """
        one_to_one_mappings = []
        for from_idx in atom_mappings:
            to_idx = atom_mappings[from_idx]
            one_to_one_mappings.append(((self.one_compound.compound_name, from_idx), (self.the_other_compound.compound_name, to_idx)))
        return one_to_one_mappings

    @staticmethod
    def detect_separate_components(compound, removed_bonds, center_atom_index):
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

        nodes = set(graph.keys() + center_atom_index)
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
    def construct_partial_compound(cpd, atom_index, removed_bonds):
        """
        To construct a partial compound based on the atom index and the removed bonds.
        :param cpd:
        :param atom_index:
        :param removed_bonds:
        :return:
        """
        atoms = [cpd.atoms[idx] for idx in atom_index]
        bonds = []
        for bond in cpd.bonds:
            atom_1, atom_2 = bond.first_atom_number, bond.second_atom_number
            if atom_1 in atom_index and atom_2 in atom_index and (atom_1, atom_2) not in removed_bonds and \
                    (atom_2, atom_1) not in removed_bonds:
                bonds.append(bond)
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
        mapped_atoms = collections.Counter()
        for atom_left in left_partial_compound.atoms:
            for atom_right in right_partial_compound.atoms:
                if atom_left.color_layers == atom_right.color_layers:
                    mapped_atoms[atom_left.atom_number] += 1
        return len(mapped_atoms) == len(left_partial_compound.atoms)

    def map_separate_components(self, left_removed_bonds, right_removed_bonds, left_index_comb, right_index_comb):
        """

        :param left_removed_bonds:
        :param right_removed_bonds:
        :param left_index_comb:
        :param right_index_comb:
        :return:
        """

        left_components = self.detect_separate_components(self.one_compound, left_removed_bonds, left_index_comb)
        right_components = self.detect_separate_components(self.the_other_compound, right_removed_bonds, right_index_comb)
        component_pairs = self.get_component_pairs(left_components, right_components)
        atom_mappings = []
        for left_component, right_component in component_pairs:
            left_partial_compound = self.construct_partial_compound(self.one_compound, left_component, left_removed_bonds)
            right_partial_compound = self.construct_partial_compound(self.the_other_compound, right_component, right_removed_bonds)
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
                collective_mappings.update(atom_mappings[idx])
            else:
                orders = groups[g]
                orders.sort(key=lambda x: self.count_different_atom_identifiers(atom_mappings[x]))
                for idx in orders:
                    if all(atom_index not in collective_mappings for atom_index in atom_mappings[idx]) and \
                            all(atom_index not in collective_mappings.values() for atom_index in atom_mappings[idx].values()):
                        collective_mappings.update(atom_mappings[idx])
        return collective_mappings

    @staticmethod
    def validate_component_atom_mappings(left_index_comb, right_index_comb, component_atom_mappings):
        """
        To check if the mapped the atoms can corresponds to the mapped reaction center atoms.
        :param left_index_comb: the list of center atom index in the left compound.
        :param right_index_comb: the list of center atom index in the right compound.
        :param component_atom_mappings: the one to one atom mappings of one separate component.
        :return:
        """

        for i, j in zip(left_index_comb, right_index_comb):
            if i in component_atom_mappings and j != component_atom_mappings[i]:
                return False
        return True

    def count_different_atom_identifiers(self, one_to_one_mappings):
        """
        To count the mapped atoms with different local atom identifier.
        :param one_to_one_mappings:
        :return:
        """
        count = 0
        for idx_1 in one_to_one_mappings:
            idx_2 = one_to_one_mappings[idx_1]
            if self.one_compound.atoms[idx_1].color_layers and self.the_other_compound.atoms[idx_2].color_layers and \
                    self.one_compound.atoms[idx_1].color_layers[1] == self.the_other_compound.atoms[idx_2].color_layers[1]:
                continue
            count += 1
        return count

def create_compounds_kcf(kcf_directory):
    """
    Create compound identities based on the KEGG kcf file.
    :param kcf_directory:
    :return:
    """
    # To put all the compounds into dictionary.
    kcf_files = glob.glob(kcf_directory + "*")
    kcf_compounds = {}
    for kcf_file in kcf_files:
        kcf_dict = kegg_kcf_parser(tools.open_text(kcf_file).split("\n"))
        atoms = [compound.Atom(atom["atom_symbol"], atom["atom_number"], x=atom["x"], y=atom["y"], kat=atom["kat"]) for
                 atom in kcf_dict["atoms"]]
        bonds = [compound.Bond(bond["first_atom_number"], bond["second_atom_number"], bond["bond_type"]) for bond in kcf_dict["bonds"]]
        kcf_compounds[kcf_dict["compound_name"]] = compound.Compound(kcf_dict["compound_name"], atoms, bonds)
    return kcf_compounds

# steps, we need to first construct KEGG kcf compounds, and use them to construct aromatic substructure set.
# then, when we construct KEGG compounds, we need to do the following steps:
# 1) check if the atom orders in the KEGG molfile and kcf are the same. and add the kat
    # we ways we can use to add the kat.
    # 1) we make use of the x, y coordinates of the atom to map between molfile and kcf.
    # 2) we use the atom coloring identifier. Here, we don't need to worry about the stereochemistry and symmetric atoms.
    # Symmetric atoms should share the same kat. In other words, atom coloring identifier is more strict than kat.
    # here I decide to use atom coloring identifier. It's safer since kegg can update the compound structure inconsistently.
# 2) add H using openbabel, the molfile should be  (for detecting bond stereochemisty)
# 2) detect aromatic substructures in the compound and update bond type.
# 3) detect the bond stereochemisty of the double bonds in the compound.
# 4) we can do the coloring or not
# 5) return the compound entity.
def create_compounds_mofile(molfile_directory, kcf_compounds, aromatic_substructures, save_file):
    """

    :param molfile_directory:
    :param kcf_compounds:
    :param aromatic_substructures:
    :return:
    """
    molfiles = glob.glob(molfile_directory + "*")
    compounds = {}
    for molfile in molfiles:
        compound_name = Path(molfile).stem
        with open(molfile, 'r') as infile:
            ct_object = ctfile.load(infile)
            # standardized the molfile representation and dump it into the openbabel to add H.
            standardized_molfile =
            # create the compound based on the ct_object.
            cpd = compound.Compound.create(ct_object, compound_name)
            # add KEGG atom type to each atom.
            cpd.color_compound()
            kcf_compound = kcf_compounds[compound_name]
            kcf_compound.color_compound()
            atom_kat = {atom.color: atom.kat for atom in kcf_compound.atoms}
            for atom in cpd.atoms:
                if atom.default_symbol != "H":
                    atom.update_kat(atom_kat[atom.color])
            # detect aromatic substructure and change aromatic bonds.
            aromatic_cycles = aromatize.detect_aromatic_substructures(cpd, aromatic_substructures)
            cpd.update_aromatic_bond_type(aromatic_cycles)
            # determine stereochemistry
            cpd.define_bond_stereochemistry()
            # return the compound.
            compounds[compound_name] = cpd
    tools.save_to_jsonpickle(compounds, save_file)
    return compounds

# when we create the kegg reaction, we need to parse the atom mappings based on rclass!
# To avoid parsing the same rclass repeatedly, let's parse the rclass first, and look it up when we need.
def create_reactions(reaction_directory, compounds, atom_mappings):

    # here we compounds, rlcass descriptions, and reactions.
    reaction_files = glob.glob(reaction_directory+"*")
    reactions = []
    for reaction_file in reaction_files:
        this_atom_mappings = []
        this_reaction = kegg_data_parser(tools.open_text(reaction_file).split("\n"))
        reaction_name = this_reaction["ENTRY"][0].split()[0]
        raw_ecs = this_reaction["ENZYME"] if "ENZYME" in this_reaction else []
        ecs = []

        for line in raw_ecs:
            ec_numbers = line.split()
            for ec in ec_numbers:
                numbers = ec.split(".")
                # some numbers in the ec are not specified. like "3.5.99.-"
                while numbers and numbers[-1] == "-":
                    numbers.pop()
                ecs.append(".".join(numbers))
        if not ecs:
            # ORTHOLOGY can contain information of ecs. Please check the example of reaction.
            for line in this_reaction["ORTHOLOGY"]:
                ec_pattern = "\[EC:.*\]"
                ec_numbers = re.findall(ec_pattern, line)[0][4:-1].split()
                for ec in ec_numbers:
                    numbers = ec.split(".")
                    while numbers and numbers[-1] == "-":
                        numbers.pop()
                    ecs.append(".".join(numbers))

        one_side_coefficients, the_other_side_coefficients = parse_equation(this_reaction["EQUATION"][0])
        one_side_compounds = [compounds[compound_name] for compound_name in one_side_coefficients]
        the_other_side_compounds = [compounds[compound_name] for compound_name in the_other_side_coefficients]
        one_side_coefficients.update(the_other_side_coefficients)
        for line in this_reaction["RCLASS"]:
            rclass, rpair = line.split()
            this_atom_mappings.extend(atom_mappings[rclass][rpair])
        reactions.append(reaction.Reaction(reaction_name, one_side_compounds, the_other_side_compounds, ecs,
                                           this_atom_mappings, one_side_coefficients))
    return reactions

def create_atom_mappings(rclass_directory, compounds):
    """

    :param rclass_directory:
    :param compounds:
    :return:
    """
    rclass_files = glob.glob(rclass_directory + "*")
    atom_mappings = collections.defaultdict(dict)
    for rclass_file in rclass_files:
        this_rclass = kegg_data_parser(tools.open_text(rclass_file).split("\n"))
        rclass_definitions = this_rclass["DEFINITION"]
        rclass_name = this_rclass["ENTRY"][0].split()[0]
        for line in this_rclass["RPAIR"]:
            tokens = line.split()
            for token in tokens:
                one_compound_name, the_other_compound_name = token.split("_")
                one_compound, the_other_compound = compounds[one_compound_name], compounds[the_other_compound_name]
                one_mappings = RpairParser(rclass_name, rclass_definitions, one_compound, the_other_compound).map_center_atoms()
                the_other_mappings = RpairParser(rclass_name, rclass_definitions, the_other_compound, one_compound).map_center_atoms()
                atom_mappings[rclass_name][token] = one_mappings if len(one_mappings) > len(the_other_mappings) else the_other_mappings
    return atom_mappings



















