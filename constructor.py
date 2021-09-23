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
    
    def __init__(self, compound_name, molfile):

        self.compound_name = compound_name
        self.molfile = molfile

    @property
    def ctfile_object(self):

        ct_object = ctfile.load(self.molfile)
        return ct_object

    def create(self):

        atoms = [ compound.Atom(atom.atom_symbol, i, atom['x'], atom['y'], atom['z'],  atom['mass_difference'], atom.charge, atom['atom_stereo_parity'],
            atom['hydrogen_count'], atom['stereo_care_box'], atom['valence'], atom['h0designator'], atom['atom_atom_mapping_number'], 
            atom['inversion_retention_flag'],  atom['exact_change_flag']) for i, atom in enumerate(self.ctfile_object.atoms) ]

        bonds = [compound.Bond(bond['first_atom_number'], bond['second_atom_number'], bond['bond_type'], bond['bond_stereo'], bond['bond_topology'],
            bond['reacting_center_status']) for bond in self.ctfile_object.bonds]

        return compound.Compound(self.compound_name, atoms, bonds)
    
        
class KCFCompoundConstructor:

    def __init__(self, kcf):

        self.kcf = kcf



class KEGGReactionConstructor:


    def __init__(self, reaction):



    def generate_atom_mappings(self, ):





class MetaCycReactionConstructor:


    def __init__(self, reaction_dict, atom_mappings):

        self.reaction_dict = reaction_dict
        self.atom_mappings = atom_mappings

    def create(self):





# class NetworkConstructor:

#     def __init__(self, compounds, reactions):

#         self.compounds = compounds
#         self.reactions = reactions
