#!/usr/bin/python3

import compound
from abc import ABC

class HarmonizedCompoundEdge:

    def __init__(self, one_compound, the_other_compound, relationship, type):

        self.one_side = one_compound
        self.the_other_side = the_other_compound
        self.relationship = relationship
        self.type = type

    # @property
    # def atom_mappings(self):
    #

class HarmonizedReactionEdge:

    def __init__(self, one_reaction, the_other_reaction, relationship):

        self.one_side = one_reaction
        self.the_other_side = the_other_reaction
        self.relationship = relationship


    # def check_atom_mappings(self):


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

        key = self.create_key(edge.one_side.name, edge.the_other_side.name)
        if key not in self.harmonized_edges:
            self.harmonized_edges[key] = edge

    def remove_edge(self, edge):

        key = self.create_key(edge.one_side.name, edge.the_other_side.name)
        if key in self.harmonized_edges:
            self.harmonized_edges.pop(key)

    def search(self, name_1, name_2):

        key = self.create_key(name_1, name_2)
        if key in self.harmonized_edges:
            return self.harmonized_edges[key]
        return None


class ReactionHarmonizationManager(HarmonizationManager):

    def __init__(self, compound_harmonization_manager):

        self.compound_harmonization_manager = compound_harmonization_manager


    @staticmethod
    def compare_ecs(one_ecs, the_other_ecs):

        if [ec for ec in one_ecs[4] if ec in the_other_ecs[4]]:
            return 4
        if [ec for ec in one_ecs[3] if ec in the_other_ecs[3]]:
            return 3
        return 0

    @staticmethod
    def harmonize_reaction(one_reaction, the_other_reaction):


def harmonize_compound_list(compound_list):

    compound_harmonization_manager = HarmonizationManager()
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
                        relationship = compounds_one[cpd_name_one].same_structure_relationship(compounds_two[cpd_name_two])
                        harmonized_compound_edge = HarmonizedCompoundEdge(compounds_one[cpd_name_one],
                                                                          compounds_two[cpd_name_two], relationship, "same_structure")
                        compound_harmonization_manager.add_edge(harmonized_compound_edge)
    return compound_harmonization_manager

def harmonize_reaction_list(reaction_list, compound_harmonization_manager):

    reaction_harmonized_manager = ReactionHarmonizationManager(compound_harmonization_manager)
    k = len(reaction_list)
    for i in range(k):
        for j in range(i+1, k):
            reactions_one, reactions_two = reaction_list[i], reaction_list[j]
            for reaction_one in reactions_one:
                for reaction_two in reactions_two:
                    reaction_harmonized_manager.harmonize_reaction(reaction_one, reaction_two)
    return reaction_harmonized_manager



