#!/usr/bin/python3
"""
mdh.harmonization
~~~~~~~~~~~~~~~~~

This module provides the :class:`~mdh.harmonization.HarmonizedEdge` class,
the :class:`~mdh.harmonization.HarmonizedCompoundEdge` class,
and the :class:`~mdh.harmonization.HarmonizedReactionEdge` class .

"""
import collections
import abc
from . import compound
from . import reaction
from typing import *


class HarmonizedEdge(abc.ABC):

    """
    The HarmonizedEdge representing compound or reaction pairs.
    """
    def __init__(self, one_side: Union[compound.Compound, reaction.Reaction],
                 the_other_side: Union[compound.Compound, reaction.Reaction], relationship: int,
                 edge_type: Union[str, int], mappings: dict) -> None:
        """
        HarmonizedEdge initializer.

        :param one_side: one side of the edge. This can be compound or reaction.
        :param the_other_side: the other side of the edge. This can be compound or reaction.
        :param relationship: equivalent, generic-specific, or loose.
        :param edge_type: for compound edge, this represents resonance, linear-circular, r group, same structure; for reaction edge, this represents 3 level match or 4 level match.
        :param mappings: for compound edge, the mappings refer to mapped atoms between compounds; for reaction edge, the mappings refer to mapped compounds between reaction.
        """
        self.one_side = one_side
        self.the_other_side = the_other_side
        self.relationship = relationship
        self.type = edge_type
        self.mappings = mappings

    @property
    def reversed_relationship(self) -> int:
        """
        To get the relationship between the other side and one side.

        :return: the reversed relationship.
        """
        if self.relationship == 2 or self.relationship == 0:
            return self.relationship
        return -self.relationship

    def pair_relationship(self, name: str) -> int:
        """
        When we map compounds in the reaction, we can access the compound edge from either side.

        :param name: the name of the searched one side.
        :return: the relationship of the searched pair.
        """
        if name == self.one_side.name:
            return self.relationship
        else:
            return self.reversed_relationship


class HarmonizedCompoundEdge(HarmonizedEdge):

    """
    The HarmonizedCompoundEdge representing compound pairs.
    """

    def __init__(self, one_compound: compound.Compound, the_other_compound: compound.Compound, relationship: int,
                 edge_type: str, atom_mappings: dict) -> None:
        """
        HarmonizedCompoundEdge initializer.

        :param one_compound: one :class:`~mdh.compound.Compound` entity in the compound pair.
        :param the_other_compound: the other :class:`~mdh.compound.Compound` entity in the compound pair.
        :param relationship: the relationship (equivalent, generic-specific, or loose) between the two compounds.
        :param edge_type: the edge type can be resonance, linear-circular, r group, or same structure.
        :param atom_mappings: the atom mappings between the two compounds.
        """
        super().__init__(one_compound, the_other_compound, relationship, edge_type, atom_mappings)

    @property
    def reversed_mappings(self) -> dict:
        """
        To get the atom mappings from compound on the other side to compound on the one side.

        :return: atom mappings between the other side compound to one side compound.
        """
        atom_mappings = collections.defaultdict(list)
        for from_atom in self.mappings:
            for to_atom in self.mappings[from_atom]:
                atom_mappings[to_atom].append(from_atom)
        return atom_mappings

    def pair_atom_mappings(self, name: str) -> dict:
        """
        To get the atom mappings of the harmonized compound edge, where one side equals to the parameter name.

        :param name: the compound name.
        :return: the atom mappings.
        """
        if name == self.one_side.name:
            return self.mappings
        else:
            return self.reversed_mappings


class HarmonizedReactionEdge(HarmonizedEdge):

    """
    The HarmonizedReactionEdge to represent reaction pairs.
    """

    def __init__(self, one_reaction: reaction.Reaction, the_other_reaction: reaction.Reaction, relationship: int,
                 edge_type: int, compound_mappings: dict) -> None:
        """
        HarmonizedReactionEdge initializer.

        :param one_reaction: one :class:`~mdh.reaction.Reaction` entity in the reaction pair.
        :param the_other_reaction: the other :class:`~mdh.reaction.Reaction` entity in the reaction pair.
        :param relationship: the relationship (equivalent, generic-specific, or loose) between the two reactions.
        :param edge_type: the reactions can be 3-level EC or 4-level EC paired.
        :param compound_mappings: the dictionary of paired compounds in the reaction pair.
        """
        super().__init__(one_reaction, the_other_reaction, relationship, edge_type, compound_mappings)

    # def check_atom_mappings(self):


class HarmonizationManager(abc.ABC):

    """
    The HarmonizationManger is responsible for adding, removing or searching harmonized edge.
    """

    def __init__(self) -> None:
        """
        HarmonizationManager initializer.
        """
        self.harmonized_edges = {}
    
    def save_manager(self) -> list:
        """
        To save all the names of harmonized edges.

        :return: the list of harmonized edges.
        """
        edges = list(self.harmonized_edges.keys())
        return edges

    @staticmethod
    def create_key(name_1: str, name_2: str) -> str:
        """
        To create the edge key. Each edge is represented by a unique key in the harmonized_edges dictionary.

        :param name_1: the name of one side of the edge.
        :param name_2: the name of the other side of the edge.
        :return: the key of the edge.
        """
        if name_1 > name_2:
            return name_1 + "@" + name_2
        else:
            return name_2 + "@" + name_1

    def add_edge(self, edge: Union[HarmonizedCompoundEdge, HarmonizedReactionEdge]) -> bool:
        """
        To add this edge to the harmonized edges.

        :param edge: the :class:`~mdh.harmonization.HarmonizedEdge` entity.
        :return: bool whether the edge is added successfully.
        """
        key = self.create_key(edge.one_side.name, edge.the_other_side.name)
        if key not in self.harmonized_edges:
            self.harmonized_edges[key] = edge
            return True
        return False

    def remove_edge(self, edge: Union[HarmonizedCompoundEdge, HarmonizedReactionEdge]) -> bool:
        """
        To remove this edge from the harmonized edges.

        :param edge: the :class:`~mdh.harmonization.HarmonizedEdge` entity.
        :return: bool whether the edge is removed successfully.
        """
        key = self.create_key(edge.one_side.name, edge.the_other_side.name)
        if key in self.harmonized_edges:
            self.harmonized_edges.pop(key)
            return True
        return False

    def search(self, name_1: str, name_2: str) -> Optional[Union[HarmonizedCompoundEdge, HarmonizedReactionEdge]]:
        """
        To search the edge based on the names of the two sides.

        :param name_1: the name of one side of the edge.
        :param name_2: the name of the other side of the edge.
        :return: edge if the edge exists or None.
        """
        key = self.create_key(name_1, name_2)
        if key in self.harmonized_edges:
            return self.harmonized_edges[key]
        return None


class CompoundHarmonizationManager(HarmonizationManager):

    """
    The CompoundHarmonizationManager is responsible for adding, removing or searching
    :class:`~mdh.harmonization.HarmonizedCompoundEdge`.
    """

    def __init__(self) -> None:
        """
        CompoundHarmonizationManager initializer.
        """
        super().__init__()
        self.compound_in_edges = collections.Counter()
        # the visited store invalid compound pairs to avoid redundant validation.
        self.visited = set()
    
    @staticmethod
    def find_compound(compound_dict: list, compound_name: str) -> Optional[compound.Compound]:
        """
        To find the :class:`~mdh.compound.Compound` based on the compound name in the compound dict.

        :param compound_dict: a list of compound dictionaries.
        :param compound_name: the target compound name.
        :return: the :class:`~mdh.compound.Compound`.
        """
        for sub_dict in compound_dict:
            if compound_name in sub_dict:
                return sub_dict[compound_name]
        return None
    
    @staticmethod
    def create_manager(compound_dict: list, compound_pairs: list):
        """
        To create the :class:`~mdh.harmonization.CompoundHarmonizationManager` based on the compound paris.

        :param compound_dict: the list of compound dictionaries.
        :param compound_pairs: the list of compound pairs.
        :return: the :class:`~mdh.harmonization.CompoundHarmonizationManager`
        """
        compound_harmonization_manager = CompoundHarmonizationManager()
        for pair in compound_pairs:
            compound_name1, compound_name2 = pair.split("@")
            compound1 = CompoundHarmonizationManager.find_compound(compound_dict, compound_name1)
            compound2 = CompoundHarmonizationManager.find_compound(compound_dict, compound_name2)
            if compound1 and compound2:
                relationship, atom_mappings = compound1.same_structure_relationship(compound2)
                harmonized_compound_edge = HarmonizedCompoundEdge(compound1, compound2, relationship, "same_structure", 
                                                                  atom_mappings)
                compound_harmonization_manager.add_edge(harmonized_compound_edge)
        return compound_harmonization_manager
                
    def add_edge(self, edge: HarmonizedCompoundEdge) -> bool:
        """
        To add a newly detected edge to the manager, and update the occurrences of compound in the harmonized edges.
        This is for calculating the jaccard index.

        :param edge: the :class:`~mdh.harmonization.HarmonizedCompoundEdge` entity.
        :return: bool whether the edge is added successfully.
        """
        if super().add_edge(edge):
            self.compound_in_edges[edge.one_side.name] += 1
            self.compound_in_edges[edge.the_other_side.name] += 1
            return True
        return False

    def remove_edge(self, edge: HarmonizedCompoundEdge) -> bool:
        """
        To remove the edge from the manager, and update the occurrences of compound in the harmonized edges.

        :param edge: the :class:`~mdh.harmonization.HarmonizedCompoundEdge` entity.
        :return: bool whether the edge is removed successfully.
        """
        if super().remove_edge(edge):
            self.compound_in_edges[edge.one_side.name] -= 1
            self.compound_in_edges[edge.the_other_side.name] -= 1
            return True
        return False

    def has_visited(self, name_1: str, name_2: str) -> bool:
        """
        To check if the compound pair has been visited.

        :param name_1: the name of one side of the edge.
        :param name_2: the name of the other side of the edge.
        :return: bool if the pair has been visited.
        """
        key = self.create_key(name_1, name_2)
        return key in self.visited

    def add_invalid(self, name_1: str, name_2: str) -> None:
        """
        To add the name of invalid compound pair to the visited.

        :param name_1: the name of one side of the edge.
        :param name_2: the name of the other side of the edge.
        :return: None.
        """
        key = self.create_key(name_1, name_2)
        self.visited.add(key)

    def get_edge_list(self) -> list:
        """
        To get the names of all the harmonized edges.

        :return: the list of names of harmonized edges.
        """
        return list(self.harmonized_edges.keys())


class ReactionHarmonizationManager(HarmonizationManager):
    """
    The ReactionHarmonizationManager is responsible for adding, removing or searching
    :class:`~mdh.harmonization.HarmonizedReactionEdge`.
    """

    def __init__(self, compound_harmonization_manager: CompoundHarmonizationManager) -> None:
        """
        ReactionHarmonizationManager initializer.

        :param compound_harmonization_manager: the :class:`~mdh.harmonization.CompoundHarmonizationManager` entity for compound pairs management.
        """
        super().__init__()
        self.compound_harmonization_manager = compound_harmonization_manager

    @staticmethod
    def compare_ecs(one_ecs: dict, the_other_ecs: dict) -> int:
        """
        To compare two lists of EC numbers.

        :param one_ecs: the dict of EC numbers of one reaction.
        :param the_other_ecs: the dict of EC numbers of the other reaction.
        :return: the level of EC number that they can be matched.
        """
        if [ec for ec in one_ecs[4] if ec in the_other_ecs[4]]:
            return 4
        if [ec for ec in one_ecs[3] if ec in the_other_ecs[3]]:
            return 3
        return 0

    @staticmethod
    def determine_relationship(relationships: list) -> int:
        """
        To determine the relationship of the reaction pair based on the relationship of paired compounds.

        :param relationships: the list of relationship of compound pairs in the two reactions.
        :return: the relationships between the two reactions.
        """
        counter = collections.Counter(relationships)
        if counter[2] > 0:
            # if at least one compound pair has loose relationship, then the relationship between the two reactions
            # is loose.
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

    def harmonize_reaction(self, one_reaction: reaction.Reaction, the_other_reaction: reaction.Reaction) -> None:
        """
        To test if two reactions can be harmonized.

        :param one_reaction: one :class:`~mdh.reaction.Reaction` that is involved in the reaction pair.
        :param the_other_reaction: the other :class:`~mdh.reaction.Reaction` that is involved in the reaction pair.
        :return: None.
        """
        ec_comparison = self.compare_ecs(one_reaction.ecs, the_other_reaction.ecs)
        if not ec_comparison:
            # Don't share the same ec, skip it.
            print("the two reactions don't share the same ec")
            return
        else:
            # For the reaction pair, one_side of one reaction can be mapped to one_side of the other reaction;
            # Or one_side of one reaction can be mapped to the_other_side of the other reaction.

            # one_side to one_side
            # we map compounds in the two sides separately.
            print("the two reactions share the same ec. Let's start comparison")
            ordered_one_side_mappings = self.compound_mappings(one_reaction.one_side, the_other_reaction.one_side)
            ordered_the_other_side_mappings = self.compound_mappings(one_reaction.the_other_side,
                                                                     the_other_reaction.the_other_side)
            ordered_jaccard = self.jaccard(one_reaction.one_side, the_other_reaction.one_side, ordered_one_side_mappings) * \
                              self.jaccard(one_reaction.the_other_side, the_other_reaction.the_other_side,
                                           ordered_the_other_side_mappings)

            # one_side to the_other_side
            reversed_one_side_mappings = self.compound_mappings(one_reaction.one_side, the_other_reaction.the_other_side)
            reversed_the_other_side_mappings = self.compound_mappings(one_reaction.the_other_side,
                                                                      the_other_reaction.one_side)
            reversed_jaccard = self.jaccard(one_reaction.one_side, the_other_reaction.the_other_side,
                                            reversed_one_side_mappings) * self.jaccard(one_reaction.the_other_side,
                                                                                       the_other_reaction.one_side,
                                                                                       reversed_the_other_side_mappings)
            # here we see which case matches better, and determine the match direction
            max_score = max(ordered_jaccard, reversed_jaccard)
            print("the jaccard score of the two reactions ", max_score)

            if ordered_jaccard > reversed_jaccard:
                one_side_pairs = [one_reaction.one_side, the_other_reaction.one_side]
                the_other_side_pairs = [one_reaction.the_other_side, the_other_reaction.the_other_side]
                mappings = [ordered_one_side_mappings, ordered_the_other_side_mappings]
            else:
                one_side_pairs = [one_reaction.one_side, the_other_reaction.the_other_side]
                the_other_side_pairs = [one_reaction.the_other_side, the_other_reaction.one_side]
                mappings = [reversed_one_side_mappings, reversed_the_other_side_mappings]

            # we need to add reaction harmonized edge, at the same time, find the compound mappings and
            # determine the relationship of the edge.
            # print(max_score)
            if max_score == 1:
                # To derive the one_to_one compound mappings on both sides.
                # print(one_reaction.name, the_other_reaction.name, mappings)
                one_side_relationships, one_side_mappings = self.one_to_one_compound_mappings(mappings[0])
                # print("one side mappings of the compound:", one_side_mappings)
                the_other_side_relationships, the_other_side_mappings = self.one_to_one_compound_mappings(mappings[1])
                # print("the other side mappings of the compound:", the_other_side_mappings)
                if one_side_relationships and the_other_side_relationships:
                    # To combine the compound mappings together.
                    one_side_mappings.update(the_other_side_mappings)
                    relationship = self.determine_relationship(one_side_relationships + the_other_side_relationships)
                    harmonized_reaction_edge = HarmonizedReactionEdge(one_reaction, the_other_reaction, relationship,
                                                                      ec_comparison, one_side_mappings)
                    print("find harmonized reaction: ", one_reaction.name, the_other_reaction.name)
                    self.add_edge(harmonized_reaction_edge)

            elif max_score >= 0.5:
                # determine if there is missed compound harmonized edge. This threshold can be adjusted.
                print("map compounds")
                one_unmapped_compounds = self.unmapped_compounds(one_side_pairs[0], one_side_pairs[1], mappings[0])
                self.match_unmapped_compounds(one_unmapped_compounds[0], one_unmapped_compounds[1])

                the_other_unmapped_compounds = self.unmapped_compounds(the_other_side_pairs[0],
                                                                       the_other_side_pairs[1], mappings[1])
                self.match_unmapped_compounds(the_other_unmapped_compounds[0], the_other_unmapped_compounds[1])

    def compound_mappings(self, one_compounds: list, the_other_compounds: list) -> dict:
        """
        To get the mapped compounds in the two compound lists.

        :param one_compounds: one list of :class:`~mdh.compound.Compound` entities.
        :param the_other_compounds: the other list of :class:`~mdh.compound.Compound` entities.
        :return: the dictionary of paired compounds with their relationship. The relationship will be used to determine the relationship of reaction pair.
        """
        mappings = collections.defaultdict(dict)
        for one_compound in one_compounds:
            for the_other_compound in the_other_compounds:
                harmonized_edge = self.compound_harmonization_manager.search(one_compound.compound_name,
                                                                             the_other_compound.compound_name)
                if harmonized_edge:
                    mappings[one_compound.compound_name][the_other_compound.compound_name] = \
                        harmonized_edge.pair_relationship(one_compound.compound_name)
        return mappings

    def unmapped_compounds(self, one_compounds: list, the_other_compounds: list, mappings: dict) -> tuple:
        """
        To get the compounds that cannot be mapped. This can lead to new compound pairs.

        :param one_compounds: one list of :class:`~mdh.compound.Compound` entities.
        :param the_other_compounds: the other list of :class:`~mdh.compound.Compound` entities.
        :param mappings: the mapped compounds between the two compound lists.
        :return: two lists of compounds that cannot be mapped.
        """
        one_side_left = [cpd for cpd in one_compounds if cpd.compound_name not in mappings.keys()]
        used_the_other = set()
        for one_cpd in mappings:
            used_the_other |= set(mappings[one_cpd].keys())
        the_other_side_left = [cpd for cpd in the_other_compounds if cpd.compound_name not in used_the_other]
        return one_side_left, the_other_side_left

    def match_unmapped_compounds(self, one_side_left: list, the_other_side_left: list) -> None:
        """
        To match the left compounds and add the valid compound pairs to the
        :class:`~mdh.harmonization.CompoundHarmonizationManager`.
        We also add the invalid compound pairs to the :class:`~mdh.harmonization.CompoundHarmonizationManager` to avoid
        redundant match.

        :param one_side_left: one list of left :class:`~mdh.compound.Compound` entities.
        :param the_other_side_left: the other list of left :class:`~mdh.compound.Compound` entities.
        :return: None.
        """
        for one_cpd in one_side_left:
            for the_other_cpd in the_other_side_left:
                print("try to map compounds, ", one_cpd.compound_name, the_other_cpd.compound_name)
                if self.compound_harmonization_manager.has_visited(one_cpd.compound_name, the_other_cpd.compound_name):
                    continue
                valid = False
                # print("current compare these two compounds: ", one_cpd.compound_name, the_other_cpd.compound_name)
                if one_cpd.formula == the_other_cpd.formula:
                    print("same formula")
                    # resonance or linear-circular type
                    try:
                        resonant_mappings = one_cpd.map_resonance(the_other_cpd, r_distance=False)
                    except:
                        self.compound_harmonization_manager.add_invalid(one_cpd.compound_name, the_other_cpd.compound_name)
                        continue
                    if resonant_mappings:
                        # here both atom_mappings and relationship should be returned.
                        print("find harmonized compound relationship with resonance relationship: ", the_other_cpd.compound_name, one_cpd.compound_name)
                        relationship, atom_mappings = one_cpd.optimal_resonant_mapping(the_other_cpd, resonant_mappings)
                        harmonized_compound_edge = HarmonizedCompoundEdge(one_cpd, the_other_cpd, relationship,
                                                                          "resonance", atom_mappings)
                        self.compound_harmonization_manager.add_edge(harmonized_compound_edge)
                        # check if they are resonant
                    else:
                        # check if they have circular and linear interchangeable formats.
                        try:
                            relationship, atom_mappings = one_cpd.circular_pair_relationship(the_other_cpd)
                        except:
                            self.compound_harmonization_manager.add_invalid(one_cpd.compound_name,
                                                                            the_other_cpd.compound_name)
                            continue
                        if atom_mappings:
                            print("find harmonized compound relationship with circular relationship: ", the_other_cpd.compound_name, one_cpd.compound_name)
                            harmonized_compound_edge = HarmonizedCompoundEdge(one_cpd, the_other_cpd, relationship,
                                                                              "circular", atom_mappings)
                            self.compound_harmonization_manager.add_edge(harmonized_compound_edge)
                        else:
                            try:
                                relationship, atom_mappings = the_other_cpd.circular_pair_relationship(one_cpd)
                            except:
                                self.compound_harmonization_manager.add_invalid(one_cpd.compound_name,
                                                                                the_other_cpd.compound_name)
                                continue
                            if atom_mappings:
                                print("find harmonized compound relationship with circular relationship: ", the_other_cpd.compound_name, one_cpd.compound_name)
                                harmonized_compound_edge = HarmonizedCompoundEdge(the_other_cpd, one_cpd, relationship,
                                                                                  "circular", atom_mappings)
                                self.compound_harmonization_manager.add_edge(harmonized_compound_edge)
                    continue

                if one_cpd.contains_r_groups():
                    # check if one cpd is more generic
                    print("try to map with r group")
                    try:
                        relationship, atom_mappings = one_cpd.with_r_pair_relationship(the_other_cpd)
                    except:
                        self.compound_harmonization_manager.add_invalid(one_cpd.compound_name,
                                                                        the_other_cpd.compound_name)
                        continue

                    if atom_mappings:
                        print("find harmonized compound relationship with r group: ", the_other_cpd.compound_name, one_cpd.compound_name)
                        harmonized_compound_edge = HarmonizedCompoundEdge(one_cpd, the_other_cpd, relationship,
                                                                          "r_group", atom_mappings)
                        self.compound_harmonization_manager.add_edge(harmonized_compound_edge)
                        valid = True

                if the_other_cpd.contains_r_groups():
                    # check if the other cpd is more generic.
                    try:
                        relationship, atom_mappings = the_other_cpd.with_r_pair_relationship(one_cpd)
                    except:
                        self.compound_harmonization_manager.add_invalid(one_cpd.compound_name,
                                                                        the_other_cpd.compound_name)
                        continue
                    if atom_mappings:
                        print("find harmonized compound relationship with r group: ", the_other_cpd.compound_name, one_cpd.compound_name)
                        harmonized_compound_edge = HarmonizedCompoundEdge(the_other_cpd, one_cpd, relationship,
                                                                          "r_group", atom_mappings)
                        self.compound_harmonization_manager.add_edge(harmonized_compound_edge)
                        valid = True
                if not valid:
                    self.compound_harmonization_manager.add_invalid(one_cpd.compound_name, the_other_cpd.compound_name)

    def jaccard(self, one_compounds: list, the_other_compounds: list, mappings: dict) -> float:
        """
        To calculate the jaccard index between the two list of compounds.

        :param one_compounds: one list of :class:`~mdh.compound.Compound` entities.
        :param the_other_compounds: the other list of :class:`~mdh.compound.Compound` entities.
        :param mappings: the dictionary of mapped compounds between the two compound lists.
        :return: the jaccard index of the two compound lists.
        """
        one_in_reactions = [self.compound_harmonization_manager.compound_in_edges[cpd.compound_name] > 0 for
                            cpd in one_compounds].count(True)
        the_other_in_reactions = [self.compound_harmonization_manager.compound_in_edges[cpd.compound_name] > 0 for cpd in
                                  the_other_compounds].count(True)
        denominator = one_in_reactions + the_other_in_reactions - len(mappings)
        if denominator > 0:
            return len(mappings) / denominator
        return 0

    def one_to_one_compound_mappings(self, mappings: dict) -> Optional[tuple]:
        """
        To find the one to one compound mappings between the two reactions.
        This step is to avoid very extreme cases that a compound in one reaction can be mapped to two or more compounds
        in the other reaction.

        :param mappings: the dictionary of compound mappings.
        :return: the tuple of relationship of compound pairs and dictionary of one to one compound mappings.
        """
        one_to_one_mappings = {}
        sorted_one_side = list(mappings.keys())
        sorted_one_side.sort()
        n = len(sorted_one_side)
        cpd_relationships = []

        def back_track(i, cpd_relationships) -> bool:
            """
            To find one to one compound mappings.

            :param i: the ith compound in the sorted_one_side.
            :param cpd_relationships: the mapped compound for the ith compound in the sorted_one_side.
            :return: bool
            """
            if i == n:
                return True
            one_side_cpd = sorted_one_side[i]
            for the_other_cpd in mappings[one_side_cpd]:
                if the_other_cpd not in one_to_one_mappings:
                    cpd_relationships.append(mappings[one_side_cpd][the_other_cpd])
                    one_to_one_mappings[the_other_cpd] = one_side_cpd
                    if not back_track(i + 1, cpd_relationships):
                        del one_to_one_mappings[the_other_cpd]
                        cpd_relationships.pop()
            return False
        back_track(0, cpd_relationships)
        # print("back tracking to achieve one to one mapping of compounds", cpd_relationships, one_to_one_mappings)
        if len(cpd_relationships) == n:
            return cpd_relationships, one_to_one_mappings
        return None, None


def harmonize_compound_list(compound_list: list) -> CompoundHarmonizationManager:
    """
    To harmonize compounds across different databases based on the compound coloring identifier.

    :param compound_list: the list of :class:`~mdh.compound.Compound` dictionary from different sources.
    :return: the :class:`~mdh.harmonization.CompoundHarmonizationManager` containing harmonized compound edges.
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
                    # here, we try to harmonize compounds with the same structure, use formula to remove unnecessary
                    # comparison.
                    # The whole compound string identifier can be rather long. Therefore, the comparison can take very
                    # long time.
                    if compounds_one[cpd_name_one].formula != compounds_two[cpd_name_two].formula:
                        continue
                    color_one = compounds_one[cpd_name_one].backbone_color_identifier(r_groups=True) + \
                                compounds_one[cpd_name_one].metal_color_identifier(details=False)
                    color_two = compounds_two[cpd_name_two].backbone_color_identifier(r_groups=True) + \
                                compounds_two[cpd_name_two].metal_color_identifier(details=False)
                    if color_one == color_two:
                        print("find one pair using color identifiers", cpd_name_one, cpd_name_two)
                        relationship, atom_mappings = compounds_one[cpd_name_one].\
                            same_structure_relationship(compounds_two[cpd_name_two])
                        harmonized_compound_edge = HarmonizedCompoundEdge(compounds_one[cpd_name_one],
                                                                          compounds_two[cpd_name_two], relationship,
                                                                          "same_structure", atom_mappings)
                        compound_harmonization_manager.add_edge(harmonized_compound_edge)
    return compound_harmonization_manager


def harmonize_reaction_list(reaction_list: list, compound_harmonization_manager: CompoundHarmonizationManager) -> \
        ReactionHarmonizationManager:
    """
    To harmonize reactions across different sources based on the harmonized compounds. At the same time, this also
    harmonizes compound pairs with resonance, linear-circular, r group types.

    :param reaction_list: a list of :class:`~mdh.reaction.Reaction` list from different sources.
    :param compound_harmonization_manager: a :class:`~mdh.harmonization.CompoundHarmonizationManager` containing harmonized compound pairs with the same structure.
    :return: :class:`~mdh.harmonization.ReactionHarmonizationManager`
    """
    reaction_harmonized_manager = ReactionHarmonizationManager(compound_harmonization_manager)
    k = len(reaction_list)
    print("k is ", k)
    last_edges = 0
    round = 0
    while len(compound_harmonization_manager.harmonized_edges) != last_edges:
        print("harmonization round {0}".format(round))
        print("current harmonized edges: ", len(compound_harmonization_manager.harmonized_edges))
        print("edges from last round ", last_edges)
        # if no new compound harmonized edges added, then we can stop harmonize.
        last_edges = len(compound_harmonization_manager.harmonized_edges)
        for i in range(k):
            for j in range(i+1, k):
                print("i and j are", i, j)
                reactions_one, reactions_two = reaction_list[i][:1000], reaction_list[j][:1000]
                print("total reaction counts ", len(reactions_one), len(reactions_two))
                for i1, reaction_one in enumerate(reactions_one):
                    for j1, reaction_two in enumerate(reactions_two):
                        print("index of two reactions, ", i1, j1, reaction_one.reaction_name, reaction_two.reaction_name)
                        reaction_harmonized_manager.harmonize_reaction(reaction_one, reaction_two)
        print("finish one round")
        round += 1
    return reaction_harmonized_manager
