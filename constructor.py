#!/usr/bin/python3

import ctfile
import compound


class CompoundConstructor:

    def __init__(mofile):

        self.molfile = molfile
    
    @property
    def ctfile(self):
 
        ctObject = ctifle.load(self.mofile)
        return ctObject
    
    def create(self):
        
        # take care of aromatic substructure detection and bond stereochemistry.

        return compound.Compound(self.ctfile)


#class KEGGCompoundConstructor(CompoundConstructor):

#    def create(self):
        
        


#class MetaCycCompoundConstructor(CompoundConstructor):

#    def create(self)



class ReactionConstructor:

    def __init__(self,)
    
    # take care of atom mapping check


class KEGGReactionConstructor:



class MetaCycReactionConstructor:



class NetworkConstructor:

    def __init__(self, compounds, reactions):

        self.compounds = compounds
        self.reactions = reactions
