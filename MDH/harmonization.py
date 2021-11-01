#!/usr/bin/python3

import collections

class HarmonizedEdge:

    def __init__(self, one_side, the_other_side, relationship, type, mappings):
        """
        
        :param one_side: 
        :param the_other_side:
        :param relationship: equivalent, generic-specific, or loose.
        :param type: for compound edge, this represents resonance, linear-circular, r group, same structure;
        for reaction edge: this represents 3 level match or 4 level match.
        """
        self.one_side = one_side
        self.the_other_side = the_other_side
        self.relationship = relationship
        self.type = type
        self.mappings = mappings

    @property
    def reversed_relationship(self):
        """

        :return:
        """
        if self.relationship == 2 or self.relationship == 0:
            return self.relationship
        return -self.relationship

    def pair_relationship(self, name):
        """

        :param name:
        :return:
        """
        if name == self.one_side.name:
            return self.relationship
        else:
            return self.reversed_relationship


class HarmonizedCompoundEdge(HarmonizedEdge):

    def __init__(self, one_compound, the_other_compound, relationship, type, atom_mappings):
        super().__init__(one_compound, the_other_compound, relationship, type, atom_mappings)

    @property
    def reversed_mappings(self):
        """
        To get the atom mappings from the_other_compound to
        :return:
        """
        atom_mappings = collections.defaultdict(list)
        for from_atom in self.mappings:
            for to_atom in self.mappings[from_atom]:
                atom_mappings[to_atom].append(from_atom)
        return atom_mappings

    def pair_atom_mappings(self, name):
        """

        :param name:
        :return:
        """
        if name == self.one_side.name:
            return self.mappings
        else:
            return self.reversed_mappings


class HarmonizedReactionEdge(HarmonizedEdge):

    def __init__(self, one_reaction, the_other_reaction, relationship, type, compound_mappings):
        super().__init__(one_reaction, the_other_reaction, relationship, type, compound_mappings)



class HarmonizationManager:

    def __init__(self):

        self.harmonized_edges = {}

    @staticmethod
    def create_key(name_1, name_2):
        """
        To create the edge key.
        :param name_1:
        :param name_2:
        :return:
        """
        if name_1 > name_2:
            return name_1 + "@" + name_2
        else:
            return name_2 + "@" + name_1

    def add_edge(self, edge):
        """

        :param edge:
        :return:
        """
        key = self.create_key(edge.one_side.name, edge.the_other_side.name)
        if key not in self.harmonized_edges:
            self.harmonized_edges[key] = edge

    def remove_edge(self, edge):
        """

        :param edge:
        :return:
        """
        key = self.create_key(edge.one_side.name, edge.the_other_side.name)
        if key in self.harmonized_edges:
            self.harmonized_edges.pop(key)

    def search(self, name_1, name_2):
        """

        :param name_1:
        :param name_2:
        :return:
        """
        key = self.create_key(name_1, name_2)
        if key in self.harmonized_edges:
            return self.harmonized_edges[key]
        return None


class CompoundHarmonizationManager(HarmonizationManager):

    def __init__(self):
        """

        """
        super().__init__()
        self.compound_in_edges = collections.Counter()

    def add_edge(self, edge):
        """

        :param edge:
        :return:
        """
        super().add_edge(edge)
        self.compound_in_edges[edge.one_side.name] += 1
        self.compound_in_edges[edge.the_other_side.name] += 1

    def remove_edge(self, edge):
        """

        :param edge:
        :return:
        """
        super().remove_edge(edge)
        self.compound_in_edges[edge.one_side.name] -= 1
        self.compound_in_edges[edge.the_other_side.name] -= 1

    def get_edge_list(self):
        """

        :return:
        """
        return list(self.harmonized_edges.keys())

class ReactionHarmonizationManager(HarmonizationManager):

    def __init__(self, compound_harmonization_manager):
        super().__init__()
        self.compound_harmonization_manager = compound_harmonization_manager

    @staticmethod
    def compare_ecs(one_ecs, the_other_ecs):
        """

        :param one_ecs:
        :param the_other_ecs:
        :return:
        """
        if [ec for ec in one_ecs[4] if ec in the_other_ecs[4]]:
            return 4
        if [ec for ec in one_ecs[3] if ec in the_other_ecs[3]]:
            return 3
        return 0

    @staticmethod
    def determine_relationship(relationships):
        """

        :param relationships:
        :return:
        """
        counter = collections.Counter(relationships)
        if counter[2] > 0:
            return 2
        if counter[1] > 0 and counter[-1] > 0:
            return 2
        if counter[1] > 0:
            return 1
        if counter[-1] > 0:
            return -1
        return 0

    def harmonize_reaction(self, one_reaction, the_other_reaction):
        """

        :param one_reaction:
        :param the_other_reaction:
        :return:
        """
        ec_comparison = self.compare_ecs(one_reaction.ecs, the_other_reaction.ecs)
        if not ec_comparison:
            # Don't share the same ec, skip it.
            return
        else:
            ordered_one_side_mappings = self.compound_mappings(one_reaction.one_side, the_other_reaction.one_side)
            ordered_the_other_side_mappings = self.compound_mappings(one_reaction.the_other_side, the_other_reaction.the_other_side)
            ordered_jaccard = self.jaccard(one_reaction.one_side, the_other_reaction.one_side, ordered_one_side_mappings) * \
                              self.jaccard(one_reaction.the_other_side, the_other_reaction.the_other_side, ordered_the_other_side_mappings)

            reversed_one_side_mappings = self.compound_mappings(one_reaction.one_side, the_other_reaction.the_other_side)
            reversed_the_other_side_mappings = self.compound_mappings(one_reaction.the_other_side, the_other_reaction.one_side)
            reversed_jaccard = self.jaccard(one_reaction.one_side, the_other_reaction.the_other_side, reversed_one_side_mappings) * \
                               self.jaccard(one_reaction.the_other_side, the_other_reaction.one_side, reversed_the_other_side_mappings)

            max_score = max(ordered_jaccard, reversed_jaccard)
            if ordered_jaccard > reversed_jaccard:
                one_side_pairs = [one_reaction.one_side, the_other_reaction.one_side]
                the_other_side_pairs = [one_reaction.the_other_side, the_other_reaction.the_other_side]
                mappings = [ordered_one_side_mappings, ordered_the_other_side_mappings]
            else:
                one_side_pairs = [one_reaction.one_side, the_other_reaction.the_other_side]
                the_other_side_pairs = [one_reaction.the_other_side, the_other_reaction.one_side]
                mappings = [reversed_one_side_mappings, reversed_the_other_side_mappings]

            # we need to add reaction harmonized edge, at the same time, find the compound mappings and determine the order.
            if  max_score == 1:
                one_side_relationships, one_side_mappings = self.compound_one_to_one_mappings(mappings[0])
                the_other_side_relationships, the_other_side_mappings = self.compound_one_to_one_mappings(mappings[1])
                if one_side_relationships and the_other_side_relationships:
                    one_side_mappings.update(the_other_side_mappings)
                    relationship = self.determine_relationship(one_side_relationships + the_other_side_relationships)
                    harmonized_reaction_edge = HarmonizedReactionEdge(one_reaction, the_other_reaction, relationship,
                                                                      ec_comparison, one_side_mappings)
                    self.add_edge(harmonized_reaction_edge)

            elif max_score >= 0.3:
            # determine if there is missed compound harmonized edge.
                one_unmapped_compounds = self.unmapped_compounds(one_side_pairs[0], one_side_pairs[1], mappings[0])
                self.match_unmapped_compounds(one_unmapped_compounds[0], one_side_pairs[1])
                the_other_unmapped_compounds = self.unmapped_compounds(the_other_side_pairs[0],
                                                                       the_other_side_pairs[1], mappings[1])
                self.match_unmapped_compounds(the_other_unmapped_compounds[0], the_other_unmapped_compounds[1])

    def unmapped_compounds(self, one_compounds, the_other_compounds, mappings):
        """

        :param one_compounds:
        :param the_other_compounds:
        :param mappings:
        :return:
        """
        one_side_left = [cpd for cpd in one_compounds if cpd.name not in mappings.keys()]
        used_the_other = set()
        for one_cpd in mappings:
            used_the_other |= set(mappings[one_cpd].keys())
        the_other_side_left = [cpd for cpd in the_other_compounds if cpd.name not in used_the_other]
        return one_side_left, the_other_side_left

    def match_unmapped_compounds(self, one_side_left, the_other_side_left):
        """

        :param one_side_left:
        :param the_other_side_left:
        :return:
        """
        for one_cpd in one_side_left:
            for the_other_cpd in the_other_side_left:
                if one_cpd.formula == the_other_cpd.formula:
                    resonant_mappings = one_cpd.map_resonance(the_other_cpd, r_distance=False)
                    if resonant_mappings:
                        # here both atom_mappings and relationship should be returned.
                        relationship, atom_mappings = one_cpd.optimal_resonant_mapping(the_other_cpd, resonant_mappings)
                        harmonized_compound_edge = HarmonizedEdge(one_cpd, the_other_cpd, relationship, "resonance", atom_mappings)
                        self.compound_harmonization_manager.add_edge(harmonized_compound_edge)
                        # check if they are resonant
                    else:
                        # check if they have circular and linear interchangeable formats.
                        relationship, atom_mappings = one_cpd.circular_pair_relationship(the_other_cpd)
                        if atom_mappings:
                            harmonized_compound_edge = HarmonizedEdge(one_cpd, the_other_cpd, relationship,"circular", atom_mappings)
                            self.compound_harmonization_manager.add_edge(harmonized_compound_edge)
                        else:
                            relationship, atom_mappings = the_other_cpd.circular_pair_relationship(one_cpd)
                            if atom_mappings:
                                harmonized_compound_edge = HarmonizedEdge(the_other_cpd, one_cpd, relationship, "circular", atom_mappings)
                                self.compound_harmonization_manager.add_edge(harmonized_compound_edge)
                    continue
                if one_cpd.contains_r_groups():
                    # check if one cpd is more generic
                    relationship, atom_mappings = one_cpd.with_r_pair_relationship(the_other_cpd)
                    if atom_mappings:
                        harmonized_compound_edge = HarmonizedEdge(the_other_cpd, one_cpd, relationship, "r_group", atom_mappings)
                        self.compound_harmonization_manager.add_edge(harmonized_compound_edge)

                if the_other_cpd.contains_r_groups():
                    # check if the other cpd is more generic.
                    relationship, atom_mappings = the_other_cpd.circular_pair_relationship(one_cpd)
                    if atom_mappings:
                        harmonized_compound_edge = HarmonizedEdge(the_other_cpd, one_cpd, relationship, "r_group", atom_mappings)
                        self.compound_harmonization_manager.add_edge(harmonized_compound_edge)

    def jaccard(self, one_compounds, the_other_compounds, mappings):
        """

        :param one_compounds:
        :param the_other_compounds:
        :param mappings:
        :return:
        """
        one_in_reactions = [self.compound_harmonization_manager.compound_in_edges[cpd.name] > 0 for cpd in one_compounds].count(True)
        the_other_in_reactions = [self.compound_harmonization_manager.compound_in_edges[cpd.name] > 0 for cpd in the_other_compounds].count(True)
        denominator = one_in_reactions + the_other_in_reactions - len(mappings)
        if denominator > 0:
            return len(mappings) / denominator
        return 0

    def compound_mappings(self, one_compounds, the_other_compounds):
        """

        :param one_compounds:
        :param the_other_compounds:
        :return:
        """
        mappings = collections.defaultdict(dict)
        for one_compound in one_compounds:
            for the_other_compound in the_other_compounds:
                harmonized_edge = self.compound_harmonization_manager.search(one_compound, the_other_compound)
                if harmonized_edge:
                    mappings[one_compound.name][the_other_compound.name] = harmonized_edge.pair_relationship(one_compound.name)
        return mappings

    def compound_one_to_one_mappings(self, mappings):
        """

        :param mappings:
        :return:
        """
        one_to_one_mappings = {}
        sorted_one_side = list(mappings.keys())
        sorted_one_side.sort()
        n = len(sorted_one_side)
        cpd_relationships = []
        def back_track(i, cpd_relationships):

            if i == n:
                return
            one_side_cpd = sorted_one_side[i]
            for the_other_cpd in mappings[one_side_cpd]:
                if the_other_cpd not in one_to_one_mappings:
                    cpd_relationships.append(mappings[one_side_cpd][the_other_cpd])
                    one_to_one_mappings[the_other_cpd] = one_side_cpd
                    back_track(i + 1, cpd_relationships)
                    del one_to_one_mappings[the_other_cpd]
                    cpd_relationships.pop()
            return
        back_track(0, cpd_relationships)
        if len(cpd_relationships) == n:
            return cpd_relationships, one_to_one_mappings
        return None, None


def harmonize_compound_list(compound_list):
    """

    :param compound_list:
    :return:
    """
    compound_harmonization_manager = CompoundHarmonizationManager()
    # here, we need to color the compound for harmonization, do it the same time to avoid redundant coloring.
    # we first just harmonize compounds with the same coloring identifiers.
    # color the compounds
    for compound_dict in compound_list:
        for compound_name in compound_dict:
            compound_dict[compound_name].color_compound(r_groups=True, bond_stereo=False, atom_stereo=False,
                                                        resonance=False, isotope_resolved=False, charge=False)
    k = len(compound_list)
    for i in range(k):
        for j in range(i + 1, k):
            compounds_one, compounds_two = compound_list[i], compound_list[j]
            for cpd_name_one in compounds_one:
                for cpd_name_two in compounds_two:
                    color_one = compounds_one[cpd_name_one].backbone_color_identifier(r_groups=True) + \
                                compounds_one[cpd_name_one].metal_color_identifier(details=False)
                    color_two = compounds_two[cpd_name_two].backbone_color_identifier(r_groups=True) + \
                                compounds_two[cpd_name_two].metal_color_identifier(details=False)
                    if color_one == color_two:
                        relationship, atom_mappings = compounds_one[cpd_name_one].same_structure_relationship(compounds_two[cpd_name_two])
                        harmonized_compound_edge = HarmonizedEdge(compounds_one[cpd_name_one],
                                                                          compounds_two[cpd_name_two], relationship,
                                                                          "same_structure", atom_mappings)
                        compound_harmonization_manager.add_edge(harmonized_compound_edge)
    return compound_harmonization_manager

def harmonize_reaction_list(reaction_list, compound_harmonization_manager):
    """

    :param reaction_list:
    :param compound_harmonization_manager:
    :return:
    """
    reaction_harmonized_manager = ReactionHarmonizationManager(compound_harmonization_manager)
    k = len(reaction_list)
    last_edges = 0

    while len(compound_harmonization_manager.harmonized_edges) != last_edges:
        last_edges = len(compound_harmonization_manager.harmonized_edges)
        for i in range(k):
            for j in range(i+1, k):
                reactions_one, reactions_two = reaction_list[i], reaction_list[j]
                for reaction_one in reactions_one:
                    for reaction_two in reactions_two:
                        reaction_harmonized_manager.harmonize_reaction(reaction_one, reaction_two)

    return reaction_harmonized_manager



