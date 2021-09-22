#!/usr/bin/python3

import ctfile
import compound
import reaction
import KEGG_parser
import MetaCyc_parser

# class CompoundConstructor:

#     def __init__(mofile):

#         self.molfile = molfile
    
#     @property
#     def ctfile(self):
 
#         ctObject = ctifle.load(self.mofile)
#         return ctObject
    
#     def create(self):
        
#         # take care of aromatic substructure detection and bond stereochemistry.

#         return compound.Compound(self.ctfile)

class MolfileCompoundConstructor:
    
    def __init__(self, compund_name, molfile):

        self.compound_name = compound_name
        self.molfile = molfile

    @property
    def ctfile(self):

        ctObject = ctifile.load(self.mofile)
        return ctObject

    def create(self):

        atoms = [ Atom(atom['x'], atom['y'], atom['z'], atom.atom_symbol, atom['mass_difference'], atom.charge, atom['atom_stereo_parity'],
            atom['hydrogen_count'], atom['stereo_care_box'], atom['valence'], atom['h0designator'], atom['atom_atom_mapping_number'], 
            atom['inversion_retention_flag'], i, atom['exact_change_flag']) for i, atom in enumeate(self.ctfile.atoms) ]

        bonds = [Bond(bond['first_atom_number'], bond['second_atom_number'], bond['bond_type'], bond['bond_stereo'], bond['bond_topology'], 
            bond['reacting_center_status']) for bond in self.ctfile.bonds]

        return compound.Compound(self.compound_name, atoms, bonds)
    
        
class KCFCompoundConstructor:

    def __init__(self, kcf):

        self.kcf = kcf

    def 




class KEGGReactionConstructor:


    def __init__(self, reaction):

        self.


    def generate_atom_mappings(self, )





class MetaCycReactionConstructor:


    def __init__(self, reaction_dict, atom_mappings):

        self.reaction_dict = reaction_dict
        self.atom_mappings = atom_mappings

    def create(self):





# class NetworkConstructor:

#     def __init__(self, compounds, reactions):

#         self.compounds = compounds
#         self.reactions = reactions
