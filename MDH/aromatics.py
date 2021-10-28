#!/usr/bin/python3

from . import compound
from . import BASS
from indigo import *

class AromaticManager:

    def __init__(self, aromatic_substructures=None):

        self.aromatic_substructures = aromatic_substructures if aromatic_substructures else []
        self.indigo = Indigo()


    def encode(self):

        return [cpd.encode() for cpd in self.aromatic_substructures]

    @staticmethod
    def decode(aromatic_structures):

        return AromaticManager([compound.Compound(sub[0], sub[1], sub[2]) for sub in aromatic_structures])

    def add_aromatic_substructures(self, substructures):

        for substructure in substructures:
            flag = False
            for aromatic_substructure in self.aromatic_substructures:
                if len(aromatic_substructure.atoms) == len(substructure.atoms) and len(aromatic_substructure.bonds) == len(substructure.bonds):
                    mapping_matrix = BASS.make_mapping_matrix(aromatic_substructure, substructure)
                    if mapping_matrix is not None:
                        isomorphs = BASS.find_mappings(aromatic_substructure.structure_matrix(resonance=False),
                                                       aromatic_substructure.distance_matrix, substructure.structure_matrix(resonance=False),
                                                       substructure.distance_matrix, mapping_matrix)
                        if isomorphs != [] and isomorphs is not None:
                            flag = True
                            break
            if not flag:
                print(substructure.compound_name, len(substructure.atoms))
                self.aromatic_substructures.append(substructure)
        return

    def kegg_aromatize(self, kcf_cpd):

        cycles = self.extract_aromatic_substructures(kcf_cpd)
        aromatic_substructures = self.construct_aromatic_entity(kcf_cpd, cycles)
        self.add_aromatic_substructures(aromatic_substructures)

    def indigo_aromatize(self, molfile):

        cpd = compound.Compound.create(molfile)
        if cpd:
            aromatic_bonds = self.indigo_aromatic_bonds(molfile)
            cycles = cpd.separate_connected_components(aromatic_bonds)
            aromatic_substructures = self.construct_aromatic_entity(cpd, cycles)
            self.add_aromatic_substructures(aromatic_substructures)

    def indigo_aromatic_bonds(self, molfile):

        aromatic_bonds = set()
        try:
            cpd = self.indigo.loadMoleculeFromFile(molfile)
            cpd.aromatize()
            for bond in cpd.iterateBonds():
                if bond.bondOrder() == 4:
                    aromatic_bonds.add((bond.source().index(), bond.destination().index()))
        except:
            print("indigo does not work for this compound")
            pass
        return aromatic_bonds

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

        aromatic_bonds = set()
        aromatic_atoms = set()
        for aromatic in self.aromatic_substructures:
            if aromatic.composition.issubset(cpd.composition):
                mapping_matrix = BASS.mappping_matrix(aromatic, cpd, True, False, False, True)
                if mapping_matrix:
                    for assignment in BASS.find_mappings(aromatic.structure_matrix(resonance=False), aromatic.distance_matrix,
                                                         cpd.structure_matrix(resonance=False), cpd.distance_matrix, mapping_matrix):
                        for bond in aromatic.bonds:
                            if aromatic.atoms[bond.first_atom_number].in_cycle and aromatic.atoms[bond.second_atom_number].in_cycle:
                                first_atom_number = cpd.heavy_atoms[assignment[bond.first_atom_number]].atom_number
                                second_atom_number = cpd.heavy_atoms[assignment[bond.second_atom_number]].atom_number
                                bond_index = (min(first_atom_number, second_atom_number), max(first_atom_number, second_atom_number))
                                aromatic_bonds.add(bond_index)
                                aromatic_atoms.add(first_atom_number)
                                aromatic_bonds.add(second_atom_number)
        cpd.update_aromatic_bond(aromatic_bonds, aromatic_atoms)

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
            bonds = [bond.clone() for bond in cpd.bonds if bond.first_atom_number in cycle and bond.second_atom_number in cycle]
            seen_atoms = set(cycle)
            for atom_index in cycle:
                atom = cpd.atoms[atom_index]
                for neighbor_index in atom.neighbors:
                    connecting_bond = cpd.bond_lookup[ (atom_index, neighbor_index)]
                    if neighbor_index not in cycle:
                        if connecting_bond.bond_type == "2":
                            bonds.append(connecting_bond.clone())
                            if neighbor_index not in seen_atoms:
                                seen_atoms.add(neighbor_index)
            atoms = [cpd.atoms[idx].clone() for idx in seen_atoms]
            idx_dict = {int(atom.atom_number): i for i, atom in enumerate(atoms)}
            for i, atom in enumerate(atoms):
                atom.update_atom_number(i)
            for bond in bonds:
                atom_1, atom_2 = bond.first_atom_number, bond.second_atom_number
                bond.update_first_atom(idx_dict[atom_1])
                bond.update_second_atom(idx_dict[atom_2])
            aromatic_substructures.append(compound.Compound(cpd.compound_name + "_" + str(count), atoms, bonds))
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

