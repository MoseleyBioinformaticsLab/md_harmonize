#!/usr/bin/python3

import compound
import BASS
from indigo import *

class AromaticManager:

    def __init__(self):

        self.aromatic_substructures = []
        self.indigo = Indigo()

    def add_aromatic_substructures(self, substructure):

        for aromatic_substructure in self.aromatic_substructures:
            if len(aromatic_substructure.atoms) == len(substructure.atoms) and len(aromatic_substructure.bonds) == len(substructure.bonds):
                mapping_matrix = BASS.make_mapping_matrix(aromatic_substructure, substructure)
                if mapping_matrix:
                    isomorphs = BASS.find_mappings(aromatic_substructure.structure_matrix, aromatic_substructure.distance_matrix,
                                                   substructure.structure_matrix, substructure.distance_matrix, mapping_matrix)
                    if isomorphs != [] and isomorphs is not None:
                        return
        self.aromatic_substructures.append(substructure)
        return

    def kegg_aromatize(self, kcf_cpd):

        cycles = self.extract_aromatic_substructures(kcf_cpd)
        aromatic_substructures = self.construct_aromatic_entity(kcf_cpd, cycles)
        for substructure in aromatic_substructures:
            self.add_aromatic_substructures(substructure)

    def indigo_aromatize(self, molfile):

        cpd = compound.Compound.create(molfile)
        aromatic_atoms = self.indigo_aromatic_atoms(molfile)
        cycles = cpd.separate_connected_components(aromatic_atoms)
        aromatic_substructures = self.construct_aromatic_entity(cpd, cycles)
        for substructure in aromatic_substructures:
            self.add_aromatic_substructures(substructure)

    def indigo_aromatic_atoms(self, molfile):

        cpd = self.indigo.loadMoleculeFromFile(molfile)
        cpd.aromatize()
        aromatic_atoms = set()
        for bond in cpd.iterateBonds():
            if bond.bondOrder() == 4:
                aromatic_atoms.add(bond.source().index())
                aromatic_atoms.add(bond.destination().index())
        return aromatic_atoms

    @staticmethod
    def fuse_cycles(cycles):
        """
        To fuse the cycles with shared atoms.
        :param cycles:
        :return:
        """
        # to remove the cycle that is contained in another cycle.
        index_set = set(index for cycle in cycles for index in cycle)
        for index in index_set:
            to_fuse = [cycle for cycle in cycles if index in cycle]
            no_fuse = [cycle for cycle in cycles if cycle not in to_fuse]
            fused = list(set([i for cycle in to_fuse for i in cycle]))
            if fused:
                cycles = no_fuse + [fused]
            if len(cycles) == 1:
                break
        return cycles

    def detect_aromatic_substructures(self, cpd):
        """
        Detect the aromatic substructures of the cpd.
        :param cpd:
        :return:
        """

        aromatic_cycles = []
        for aromatic in self.aromatic_substructures:
            if aromatic.composition.issubset(cpd.composition):
                mmat = BASS.mappping_matrix(aromatic, cpd, True, False, False, True)
                if mmat:
                    for
        aromatic_cycles = self.fuse_cycles(aromatic_cycles)
        cpd.update_aromatic_bond_type(aromatic_cycles)

    @staticmethod
    def construct_aromatic_entity(cpd, aromatic_cycles):
        """
        To construct the aromatic substructure compound entity based on the aromatic atoms.
        :param cpd:
        :param aromatic_cycles:
        :return:
        """
        count = 0
        aromatic_substructures = []
        for cycle in aromatic_cycles:
            bonds = [bond for bond in cpd.bonds if bond.first_atom_number in cycle and bond.second_atom_number in cycle]
            atoms = []
            seen_atoms = set(cycle)
            for atom_index in cycle:
                atom = cpd.atoms[atom_index]
                atoms.append(atom)
                for neighbor_index in atom.neighbors:
                    connecting_bond = cpd.bond_lookup[ (atom_index, neighbor_index)]
                    if neighbor_index not in cycle:
                        if connecting_bond.bond_type == "2":
                            bonds.append(connecting_bond)
                            if neighbor_index not in seen_atoms:
                                seen_atoms.add(neighbor_index)
                                atoms.append(cpd.atoms[neighbor_index])
                aromatic_substructures.append(compound.Compound(cpd.compound_name + str(count), atoms, bonds))
            count += 1
        return aromatic_substructures

    def extract_aromatic_substructures(self, cpd):
        """
        To detect the aromatic substructures in a compound based on the aromatic types. This is just for KEGG kcf file.
        :param cpd:
        :return:
        """
        aromatic_elements = ["C", "N"]
        aromatic_types = ["C8x", "C8y", "N4x", "N4y", "N5x", "N5y"]
        aromatic_cycles = []
        for cycle in cpd.cycles:
            aromatic_atoms = [cpd.atoms[index] for index in cycle if cpd.atoms[index].default_symbol in aromatic_elements]
            is_aromatic = [ atom.kat in aromatic_types for atom in aromatic_atoms ].count(True) == len(aromatic_atoms)
            oxygen_count = [cpd.atoms[index].default_symbol for index in cycle].count("O")
            # see the middle in KEGG compound C03861.
            if not is_aromatic or (oxygen_count >= 2 and len(cycle) >= 6):
                continue
            aromatic_cycles.append(cycle)
        return self.fuse_cycles(aromatic_cycles)








