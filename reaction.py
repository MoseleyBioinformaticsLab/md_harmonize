#!/usr/bin/python3

import json


class Edge:


    def __init__(self, one_compound, one_atom_idx, the_other_compound, the_other_atom_idx):
        """

        """
        self.one_compound = one_compound
        self.one_atom_idx = one_atom_idx
        self.the_other_compound = the_other_compound
        self.the_other_atom_idx


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
        """
    
        """
        