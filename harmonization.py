#!/usr/bin/python3
import compound
from abc import ABC

class HarmonizedCompoundEdge:

    def __init__(self, one_compound, the_other_compound, relationship):

        self.one_side = one_compound
        self.the_other_side = the_other_compound
        self.relationship = relationship

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





