#!/usr/bin/python3
"""
MDH.harmonization
~~~~~~~~~~~~~~~~~

This module provides the :class:`~MDH.harmonization.HarmonizedEdge` class, the :class:`MDH.harmonization.HarmonizedCompoundEdge` class,
and the :class:`~MDH.harmonization.HarmonizedReactionEdge` class .

"""
import collections
import abc

class HarmonizedEdge(abc.ABC):

    """
    The HarmonizedEdge to represent compound or reaction pairs.
    """

    def __init__(self, one_side, the_other_side, relationship, edge_type, mappings):
        """
        HarmonizedEdge initializer.

        :param one_side: one side of the edge. This can be compound or reaction.
        :type one_side: :class:`~MDH.compound.Compound` or :class:`~MDH.reaction.Reaction`.
        :param the_other_side: the other side of the edge. This can be compound or reaction.
        :type the_other_side: :class:`~MDH.compound.Compound` or :class:`~MDH.reaction.Reaction`.
        :param relationship: equivalent, generic-specific, or loose.
        :type relationship: :py:class:`int`.
        :param edge_type: for compound edge, this represents resonance, linear-circular, r group, same structure;
        for reaction edge, this represents 3 level match or 4 level match.
        :type edge_type: :py:class:`str` or :py:class:`int`.
        :param mappings: for compound edge, the mappings refer to mapped atoms between compounds; for reaction edge,
        the mappings refer to mapped compounds between reaction.
        :type mappings: :py:class:`dict`.
        """
        self.one_side = one_side
        self.the_other_side = the_other_side
        self.relationship = relationship
        self.type = edge_type
        self.mappings = mappings

    @property
    def reversed_relationship(self):
        """
        To get the relationship between the other side and one side.

        :return: the reversed relationship.
        :rtype: :py:class:`int`.
        """
        if self.relationship == 2 or self.relationship == 0:
            return self.relationship
        return -self.relationship

    def pair_relationship(self, name):
        """
        When we map compounds in the reaction, we can access the compound edge from either side.

        :param name: the name of the searched one side.
        :type name: :py:class:`str`.
        :return: the relationship of the searched pair.
        :rtype: :py:class:`int`.
        """
        if name == self.one_side.name:
            return self.relationship
        else:
            return self.reversed_relationship


class HarmonizedCompoundEdge(HarmonizedEdge):

    """
    The HarmonizedCompoundEdge to represent compound pairs.
    """

    def __init__(self, one_compound, the_other_compound, relationship, edge_type, atom_mappings):
        """
        HarmonizedCompoundEdge initializer.

        :param one_compound: one :class:`~MDH.compound.Compound` entity in the compound pair.
        :type one_compound: :class:`~MDH.compound.Compound`
        :param the_other_compound: the other :class:`~MDH.compound.Compound` entity in the compound pair.
        :type the_other_compound: :class:`~MDH.compound.Compound`
        :param relationship: the relationship (equivalent, generic-specific, or loose) between the two compounds.
        :type relationship: :py:class:`int`.
        :param edge_type: the edge type can be resonance, linear-circular, r group, or same structure.
        :type edge_type: :py:class:`str`
        :param atom_mappings: the atom mappings between the two compounds.
        :type atom_mappings: :py:class:`dict`.
        """
        super().__init__(one_compound, the_other_compound, relationship, edge_type, atom_mappings)

    @property
    def reversed_mappings(self):
        """
        To get the atom mappings from compound on the other side to compound on the one side.

        :return: atom mappings between the other side compound to one side compound.
        :rtype: :py:class:`dict`.
        """
        atom_mappings = collections.defaultdict(list)
        for from_atom in self.mappings:
            for to_atom in self.mappings[from_atom]:
                atom_mappings[to_atom].append(from_atom)
        return atom_mappings

    def pair_atom_mappings(self, name):
        """
        To get the atom mappings of the harmonized compound edge, where one side equals to the parameter name.

        :param name: the compound name.
        :type name: :py:class:`str`.
        :return: the atom mappings.
        :rtype: :py:class:`dict`.
        """
        if name == self.one_side.name:
            return self.mappings
        else:
            return self.reversed_mappings


class HarmonizedReactionEdge(HarmonizedEdge):

    """
    The HarmonizedReactionEdge to represent reaction pairs.
    """

    def __init__(self, one_reaction, the_other_reaction, relationship, edge_type, compound_mappings):
        """
        HarmonizedReactionEdge initializer.

        :param one_reaction: one :class:`~MDH.reaction.Reaction` entity in the reaction pair.
        :type one_reaction: :class:`~MDH.reaction.Reaction`.
        :param the_other_reaction: the other :class:`~MDH.reaction.Reaction` entity in the reaction pair.
        :type the_other_reaction: :class:`~MDH.reaction.Reaction`.
        :param relationship: the relationship (equivalent, generic-specific, or loose) between the two reactions.
        :type relationship: :py:class:`str`.
        :param edge_type: the reactions can be 3-level EC or 4-level EC paired.
        :type edge_type: :py:class:`int`.
        :param compound_mappings: the dictionary of paired compounds in the reaction pair.
        :rtype compound_mappings: :py:class:`dict`.
        """
        super().__init__(one_reaction, the_other_reaction, relationship, edge_type, compound_mappings)

    # def check_atom_mappings(self):


class HarmonizationManager(abc.ABC):

    """
    The HarmonizationManger takes charge of adding, removing or searching harmonized edge.
    """

    def __init__(self):
        """
        HarmonizationManager initializer.
        """
        self.harmonized_edges = {}

    @staticmethod
    def create_key(name_1, name_2):
        """
        To create the edge key. Each edge is represented by a unique key in the harmonized_edges dictionary.

        :param name_1: the name of one side of the edge.
        :type name_1: :py:class:`str`.
        :param name_2: the name of the other side of the edge.
        :type name_2: :py:class:`str`.
        :return: the key of the edge.
        :rtype: :py:class:`str`.
        """
        if name_1 > name_2:
            return name_1 + "@" + name_2
        else:
            return name_2 + "@" + name_1

    def add_edge(self, edge):
        """
        To add this edge to the harmonized edges.

        :param edge: the :class:`~MDH.harmonization.HarmonizedEdge` entity.
        :type edge: :class:`~MDH.harmonization.HarmonizedEdge`.
        :return: bool whether the edge does not exist and is added successfully.
        :rtype: :py:obj:`bool`.
        """
        key = self.create_key(edge.one_side.name, edge.the_other_side.name)
        if key not in self.harmonized_edges:
            self.harmonized_edges[key] = edge
            return True
        return False

    def remove_edge(self, edge):
        """
        To remove this edge from the harmonized edges.

        :param edge: the :class:`~MDH.harmonization.HarmonizedEdge` entity.
        :type edge: :class:`~MDH.harmonization.HarmonizedEdge`.
        :return: bool whether the edge exists and is removed successfully.
        :rtype: :py:obj:`bool`.
        """
        key = self.create_key(edge.one_side.name, edge.the_other_side.name)
        if key in self.harmonized_edges:
            self.harmonized_edges.pop(key)
            return True
        return False

    def search(self, name_1, name_2):
        """
        To search the edge based on the names of the two sides.

        :param name_1: the name of one side of the edge.
        :type name_1: :py:class:`str`.
        :param name_2: the name of the other side of the edge.
        :type name_2: :py:class:`str`.
        :return: edge if the edge exists or None.
        :rtype: :class:`~MDH.harmonization.HarmonizedEdge` or :py:obj:`None`.
        """
        key = self.create_key(name_1, name_2)
        if key in self.harmonized_edges:
            return self.harmonized_edges[key]
        return None


class CompoundHarmonizationManager(HarmonizationManager):

    """
    The CompoundHarmonizationManager takes charge of adding, removing or searching :class:`~MDH.harmonization.HarmonizedCompoundEdge`.
    """

    def __init__(self):
        """
        CompoundHarmonizationManager initializer.
        """
        super().__init__()
        self.compound_in_edges = collections.Counter()
        # the visited store invalid compound pairs to avoid redundant validation.
        self.visited = set()

    def add_edge(self, edge):
        """
        To add a newly detected edge to the manager, and update the occurrences of compound in the harmonized edges.
        This is for calculating the jaccard index.

        :param edge: the :class:`~MDH.harmonization.HarmonizedCompoundEdge` entity.
        :type edge: :class:`~MDH.harmonization.HarmonizedCompoundEdge`.
        :return: bool whether the edge does not exist and is added successfully.
        :rtype: :py:obj:`bool`.
        """
        if super().add_edge(edge):
            self.compound_in_edges[edge.one_side.name] += 1
            self.compound_in_edges[edge.the_other_side.name] += 1
            return True
        return False

    def remove_edge(self, edge):
        """
        To remove the edge from the manager, and update the occurrences of compound in the harmonized edges.

        :param edge: the :class:`~MDH.harmonization.HarmonizedCompoundEdge` entity.
        :type edge: :class:`~MDH.harmonization.HarmonizedCompoundEdge`.
        :return: bool whether the edge exists and is removed successfully.
        :rtype: :py:obj:`bool`.
        """
        if super().remove_edge(edge):
            self.compound_in_edges[edge.one_side.name] -= 1
            self.compound_in_edges[edge.the_other_side.name] -= 1
            return True
        return False

    def has_visited(self, name_1, name_2):
        """
        To check if the compound pair has been visited.

        :param name_1: the name of one side of the edge.
        :type name_1: :py:class:`str`.
        :param name_2: the name of the other side of the edge.
        :type name_2: :py:class:`str`.
        :return: bool if the pair has been visited.
        :rtype: :py:obj:`bool`.
        """
        key = self.create_key(name_1, name_2)
        return key in self.visited

    def add_invalid(self, name_1, name_2):
        """
        To add the name of invalid compound pair to the visited.

        :param name_1: the name of one side of the edge.
        :type name_1: :py:class:`str`.
        :param name_2: the name of the other side of the edge.
        :type name_2: :py:class:`str`.
        :return: None.
        :rtype: :py:obj:`None`.
        """
        key = self.create_key(name_1, name_2)
        self.visited.add(key)

    def get_edge_list(self):
        """
        To get the names of all the harmonized edges.

        :return: the list of names of harmonized edges.
        :rtype: :py:class:`list`.
        """
        return list(self.harmonized_edges.keys())

class ReactionHarmonizationManager(HarmonizationManager):
    """
    The ReactionHarmonizationManager takes charge of adding, removing or searching :class:`~MDH.harmonization.HarmonizedReactionEdge`.
    """

    def __init__(self, compound_harmonization_manager):
        """
        ReactionHarmonizationManager initializer.

        :param compound_harmonization_manager: the :class:`~MDH.harmonization.CompoundHarmonizationManager` entity for
        compound pairs management.
        :type compound_harmonization_manager: :class:`~MDH.harmonization.CompoundHarmonizationManager`.
        """
        super().__init__()
        self.compound_harmonization_manager = compound_harmonization_manager

    @staticmethod
    def compare_ecs(one_ecs, the_other_ecs):
        """
        To compare two lists of EC numbers.

        :param one_ecs: one list of EC numbers.
        :type one_ecs: :py:class:`list`.
        :param the_other_ecs: the other list of EC numbers.
        :type the_other_ecs: :py:class:`list`.
        :return: the level of EC number that they can matched.
        :rtype: :py:class:`int`.
        """
        if [ec for ec in one_ecs[4] if ec in the_other_ecs[4]]:
            return 4
        if [ec for ec in one_ecs[3] if ec in the_other_ecs[3]]:
            return 3
        return 0

    @staticmethod
    def determine_relationship(relationships):
        """
        To determine the relationship of the reaction pair based on the relationship of paired compounds.

        :param relationships: the list of relationship of compound pairs in the two reactions.
        :type relationships: :py:class:`list`.
        :return: the relationships between the two reactions.
        :rtype: :py:class:`int`.
        """
        counter = collections.Counter(relationships)
        if counter[2] > 0:
            # if at least one compound pair has loose relationship, then the relationship between the two reactions is loose.
            return 2
        if counter[1] > 0 and counter[-1] > 0:
            # if compound pairs can have generic-specific as well specific-generic relationship, then the relationship
            # between the two reactions is loose.
            return 2
        if counter[1] > 0:
            return 1
        if counter[-1] > 0:
            return -1
        # equivalent relationship.
        return 0

    def harmonize_reaction(self, one_reaction, the_other_reaction):
        """
        To test if two reactions can be harmonized.

        :param one_reaction: one :class:`~MDH.reaction.Reaction` that is involved in the reaction pair.
        :type one_reaction: :class:`~MDH.reaction.Reaction`.
        :param the_other_reaction: the other :class:`~MDH.reaction.Reaction` that is involved in the reaction pair.
        :type the_other_reaction: :class:`~MDH.reaction.Reaction`.
        :return: None.
        :rtype: :py:obj:`None`.
        """
        ec_comparison = self.compare_ecs(one_reaction.ecs, the_other_reaction.ecs)
        if not ec_comparison:
            # Don't share the same ec, skip it.
            return
        else:
            # For the reaction pair, one_side of one reaction can be mapped to one_side of the other reaction;
            # Or one_side of one reaction can be mapped to the_other_side of the other reaction.

            # one_side to one_side
            # we map compounds in the two sides separately.
            ordered_one_side_mappings = self.compound_mappings(one_reaction.one_side, the_other_reaction.one_side)
            ordered_the_other_side_mappings = self.compound_mappings(one_reaction.the_other_side, the_other_reaction.the_other_side)
            ordered_jaccard = self.jaccard(one_reaction.one_side, the_other_reaction.one_side, ordered_one_side_mappings) * \
                              self.jaccard(one_reaction.the_other_side, the_other_reaction.the_other_side, ordered_the_other_side_mappings)

            #one_side to the_other_side
            reversed_one_side_mappings = self.compound_mappings(one_reaction.one_side, the_other_reaction.the_other_side)
            reversed_the_other_side_mappings = self.compound_mappings(one_reaction.the_other_side, the_other_reaction.one_side)
            reversed_jaccard = self.jaccard(one_reaction.one_side, the_other_reaction.the_other_side, reversed_one_side_mappings) * \
                               self.jaccard(one_reaction.the_other_side, the_other_reaction.one_side, reversed_the_other_side_mappings)


            # here we see which case matches better, and determine the match direction
            max_score = max(ordered_jaccard, reversed_jaccard)
            if ordered_jaccard > reversed_jaccard:
                one_side_pairs = [one_reaction.one_side, the_other_reaction.one_side]
                the_other_side_pairs = [one_reaction.the_other_side, the_other_reaction.the_other_side]
                mappings = [ordered_one_side_mappings, ordered_the_other_side_mappings]
            else:
                one_side_pairs = [one_reaction.one_side, the_other_reaction.the_other_side]
                the_other_side_pairs = [one_reaction.the_other_side, the_other_reaction.one_side]
                mappings = [reversed_one_side_mappings, reversed_the_other_side_mappings]

            # we need to add reaction harmonized edge, at the same time, find the compound mappings and determine the relationship of the edge.
            if  max_score == 1:
                # To derive the one_to_one compound mappings on both sides.
                one_side_relationships, one_side_mappings = self.one_to_one_compound_mappings(mappings[0])
                the_other_side_relationships, the_other_side_mappings = self.one_to_one_compound_mappings(mappings[1])
                if one_side_relationships and the_other_side_relationships:
                    # To combine the compound mappings together.
                    one_side_mappings.update(the_other_side_mappings)
                    relationship = self.determine_relationship(one_side_relationships + the_other_side_relationships)
                    harmonized_reaction_edge = HarmonizedReactionEdge(one_reaction, the_other_reaction, relationship,
                                                                      ec_comparison, one_side_mappings)
                    self.add_edge(harmonized_reaction_edge)

            elif max_score >= 0.3:
            # determine if there is missed compound harmonized edge. This threshold can be adjusted.
                one_unmapped_compounds = self.unmapped_compounds(one_side_pairs[0], one_side_pairs[1], mappings[0])
                self.match_unmapped_compounds(one_unmapped_compounds[0], one_side_pairs[1])

                the_other_unmapped_compounds = self.unmapped_compounds(the_other_side_pairs[0],
                                                                       the_other_side_pairs[1], mappings[1])
                self.match_unmapped_compounds(the_other_unmapped_compounds[0], the_other_unmapped_compounds[1])

    def compound_mappings(self, one_compounds, the_other_compounds):
        """
        To get the mapped compounds in the two compound lists.

        :param one_compounds: one list of :class:`~MDH.compound.Compound` entities.
        :type one_compounds: :py:class:`list`.
        :param the_other_compounds: the other list of :class:`~MDH.compound.Compound` entities.
        :type the_other_compounds: :py:class:`list`.
        :return: the dictionary of paired compounds with their relationship. The relationship will be used to determine
        the relationship of reaction pair.
        :rtype: :py:class:`dict`.
        """
        mappings = collections.defaultdict(dict)
        for one_compound in one_compounds:
            for the_other_compound in the_other_compounds:
                harmonized_edge = self.compound_harmonization_manager.search(one_compound, the_other_compound)
                if harmonized_edge:
                    mappings[one_compound.name][the_other_compound.name] = harmonized_edge.pair_relationship(one_compound.name)
        return mappings

    def unmapped_compounds(self, one_compounds, the_other_compounds, mappings):
        """
        To get the compounds that cannot be mapped. This can lead to new compound pairs.

        :param one_compounds: one list of :class:`~MDH.compound.Compound` entities.
        :type one_compounds: :py:class:`list`.
        :param the_other_compounds: the other list of :class:`~MDH.compound.Compound` entities.
        :type the_other_compounds: :py:class:`list`.
        :param mappings: the mapped compounds between the two compound lists.
        :type mappings: :py:class:`dict`.
        :return: two lists of compounds that cannot be mapped.
        :rtype: :py:class:`list`.
        """
        one_side_left = [cpd for cpd in one_compounds if cpd.name not in mappings.keys()]
        used_the_other = set()
        for one_cpd in mappings:
            used_the_other |= set(mappings[one_cpd].keys())
        the_other_side_left = [cpd for cpd in the_other_compounds if cpd.name not in used_the_other]
        return one_side_left, the_other_side_left

    def match_unmapped_compounds(self, one_side_left, the_other_side_left):
        """
        To match the left compounds and add the valid compound pairs to the :class:`~MDH.harmonization.CompoundHarmonizationManager`.
        We also add the invalid compound pairs to the :class:`~MDH.harmonization.CompoundHarmonizationManager` to avoid
        redundant match.

        :param one_side_left: one list of left :class:`~MDH.compound.Compound` entities.
        :type one_side_left: :py:class:`list`.
        :param the_other_side_left: the other list of left :class:`~MDH.compound.Compound` entities.
        :type the_other_side_left: :py:class:`list`.
        :return: None.
        :rtype: :py:obj:`None`.
        """
        for one_cpd in one_side_left:
            for the_other_cpd in the_other_side_left:
                if self.compound_harmonization_manager.has_visited(one_cpd.name, the_other_cpd.name):
                    continue
                valid = False
                if one_cpd.formula == the_other_cpd.formula:
                    # resonance or linear-circular type
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
                        valid = True

                if the_other_cpd.contains_r_groups():
                    # check if the other cpd is more generic.
                    relationship, atom_mappings = the_other_cpd.circular_pair_relationship(one_cpd)
                    if atom_mappings:
                        harmonized_compound_edge = HarmonizedEdge(the_other_cpd, one_cpd, relationship, "r_group", atom_mappings)
                        self.compound_harmonization_manager.add_edge(harmonized_compound_edge)
                        valid = True

                if not valid:
                    self.compound_harmonization_manager.add_invalid(one_cpd.name, the_other_cpd.name)

    def jaccard(self, one_compounds, the_other_compounds, mappings):
        """
        To calculate the jaccard index between the two list of compounds.

        :param one_compounds: one list of :class:`~MDH.compound.Compound` entities.
        :type one_compounds: :py:class:`list`.
        :param the_other_compounds: the other list of :class:`~MDH.compound.Compound` entities.
        :type the_other_compounds: :py:class:`list`.
        :param mappings: the dictionary of mapped compounds between the two compound lists.
        :type mappings: :py:class:`dict`.
        :return: the jaccard index of the two compound lists.
        :rtype: :py:class:`float`.
        """
        one_in_reactions = [self.compound_harmonization_manager.compound_in_edges[cpd.name] > 0 for cpd in one_compounds].count(True)
        the_other_in_reactions = [self.compound_harmonization_manager.compound_in_edges[cpd.name] > 0 for cpd in the_other_compounds].count(True)
        denominator = one_in_reactions + the_other_in_reactions - len(mappings)
        if denominator > 0:
            return len(mappings) / denominator
        return 0

    def one_to_one_compound_mappings(self, mappings):
        """
        To find the one to one compound mappings between the two reactions.
        This step is to avoid very extreme cases that a compound in one reaction can be mapped to two or more compounds in the other reaction.

        :param mappings: the dictionary of compound mappings.
        :type mappings: :py:class:`dict`.
        :return: the list of relationship of compound pairs and dictionary of one to one compound mappings.
        :rtype: :py:class:`list` and :py:class:`dict`..
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
    To harmonize compounds across different databases based on the compound coloring identifier.

    :param compound_list: the list of :class:`~MDH.compound.Compound` dictionary from different sources.
    :type compound_list: :py:class:`list`.
    :return: the :class:`~MDH.harmonization.CompoundHarmonizationManager` entity with harmonized compound edges.
    :rtype: :class:`~MDH.harmonization.CompoundHarmonizationManager`.
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
    To harmonize reactions across different sources based on the harmonized compounds. At the same time, this also harmonizes
    compound pairs with resonance, linear-circular, r group types.

    :param reaction_list: a list of :class:`~MDH.reaction.Reaction` list from different sources.
    :type reaction_list: :py:class:`list`.
    :param compound_harmonization_manager: a :class:`~MDH.harmonization.CompoundHarmonizationManager` containing
    harmonized compound pairs with the same structure.
    :return: :class:`~MDH.harmonization.ReactionHarmonizationManager`
    """
    reaction_harmonized_manager = ReactionHarmonizationManager(compound_harmonization_manager)
    k = len(reaction_list)
    last_edges = 0

    while len(compound_harmonization_manager.harmonized_edges) != last_edges:
        # if no new compound harmonized edges added, then we can stop harmonize.
        last_edges = len(compound_harmonization_manager.harmonized_edges)
        for i in range(k):
            for j in range(i+1, k):
                reactions_one, reactions_two = reaction_list[i], reaction_list[j]
                for reaction_one in reactions_one:
                    for reaction_two in reactions_two:
                        reaction_harmonized_manager.harmonize_reaction(reaction_one, reaction_two)

    return reaction_harmonized_manager



