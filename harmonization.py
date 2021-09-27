#!/usr/bin/python3
import compound

class HarmonizedCompoundEdge:

    def __init__(self, one_compound, the_other_compound, relationship):

        self.one_compound = one_compound
        self.the_other_compound = the_other_compound
        self.relationship = relationship

    # @property
    # def atom_mappings(self):
    #

class HarmonizedReactionEdge:

    def __init__(self, one_reaction, the_other_reaction, relationship):

        self.one_reaction = one_reaction
        self.the_other_reaction = the_other_reaction
        self.relationship = relationship


    # def check_atom_mappings(self):


def compound_harmonization(compound_list_1, compound_list_2, save_file):
    """

    :param compound_list_1:
    :param compound_list_2:
    :return:
    """
    harmonized_compound_edges = []
    for cpd_1 in compound_list_1:
        for cpd_2 in compound_list_2:


    return harmonized_compound_edges

def reaction_harmonization(reaction_list_1, reaction_list_2, save_file):
    """

    :param reaction_list_1:
    :param reaction_list_2:
    :return:
    """
    harmonized_reaction_edges = []
    for reaction_1 in reaction_list_1:
        for reaction_2 in reaction_list_2:

    return harmonized_reaction_edges

