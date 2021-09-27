#!/usr/bin/python3

import json


class ReactionCompoundEdge:


    def __init__(self, one_compound, one_atom_index, the_other_compound, the_other_atom_index):
        """

        """
        self.one_compound = one_compound
        self.one_atom_index = one_atom_index
        self.the_other_compound = the_other_compound
        self.the_other_atom_index = the_other_atom_index


class Reaction:

    def __init__(self, reaction_name, one_compound, the_other_compound, ECs, atom_mappings):
        
        """
        

        """
        self.reaction_name = reaction_name
        self.one_compound = one_compound
        self.the_other_compound = the_other_compound
        self.ECs = ECs 
        self.atom_mappings = atom_mappings
        
    
    def create_edge(self):

    
    def find_reaction_pair(self, the_other):

        