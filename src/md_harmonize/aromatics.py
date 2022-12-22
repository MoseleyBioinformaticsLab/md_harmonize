#!/usr/bin/python3

"""
md_harmonize.aromatics
~~~~~~~~~~~~~~~~~~~~~~

This module provides the :class:`~md_harmonize.aromatics.AromaticManager` class entity.

"""

from . import tools
from . import compound
from indigo import *
try:
    from . import BASS
except ImportError:
    from . import uncythonized_BASS as BASS

class AromaticManager:

    """
    Two major functions are implemented in AromaticManager.

    1) Extract aromatic substructures based on labelled aromatic atoms (mainly C and N) or Indigo detected aromatic bonds.
        The first case only applies to KEGG compounds, and the second case applies to compounds from any databases.

    2) Detect the aromatic substructures in any given compound, and update the bond type of the detected aromatic bonds.
    """

    def __init__(self, aromatic_substructures: list = None) -> None:
        """
        AromaticManager initializer.

        :param aromatic_substructures: a list of aromatic substructures.
        """
        self.aromatic_substructures = aromatic_substructures if aromatic_substructures else []
        self.indigo = Indigo()

        for substructure in self.aromatic_substructures:
            substructure.update_color_tuple()

    def encode(self) -> list:
        """
        To encode the aromatic substructures in the aromatic manager. (Get error when try to jsonpickle the
        AromaticManager: the cythonized entities cannot be pickled.)

        :return: the list of aromatic substructures.
        """
        return [cpd.encode() for cpd in self.aromatic_substructures]

    @staticmethod
    def decode(aromatic_structures: list):
        """
        Construct the AromaticManager based on the aromatic substructures.

        :param aromatic_structures: the list of aromatic substructures.
        :return: the constructed AromaticManager.
        """
        return AromaticManager([compound.Compound(sub[0], sub[1], sub[2]) for sub in aromatic_structures])

    def add_aromatic_substructures(self, substructures: list) -> None:
        """
        Add newly detected aromatic structures to the manager. Make sure no duplicates in the aromatic substructures.

        :param substructures: a list of aromatic substructures.
        :return: None.
        """
        for substructure in substructures:
            flag = False
            for aromatic_substructure in self.aromatic_substructures:
                if all(aromatic_substructure.composition[key] == substructure.composition[key] for key in
                       substructure.composition):
                    mapping_matrix = BASS.make_mapping_matrix(aromatic_substructure, substructure, True, True, False)
                    if mapping_matrix is not None:
                        isomorphs = BASS.find_mappings(aromatic_substructure.structure_matrix(resonance=False),
                                                       aromatic_substructure.distance_matrix,
                                                       substructure.structure_matrix(resonance=False),
                                                       substructure.distance_matrix, mapping_matrix)
                        if isomorphs != [] and isomorphs is not None:
                            flag = True
                            break
            if not flag:
                substructure.update_color_tuple()
                self.aromatic_substructures.append(substructure)
        return

    def kegg_aromatize(self, kcf_cpd: compound.Compound) -> None:
        """
        Extract aromatic substructures based on KEGG atom type in KEGG compound parsed from KCF file, and add the
        newly detected aromatic substructures to the AromaticManager.

        :param kcf_cpd: the KEGG compound entity derived from KCF file.
        :return: None.
        """
        cycles = self.extract_aromatic_substructures(kcf_cpd)
        aromatic_substructures = self.construct_aromatic_entity(kcf_cpd, cycles)
        self.add_aromatic_substructures(aromatic_substructures)

    def indigo_aromatize(self, molfile: str) -> None:
        """
        Extract aromatic substructures via Indigo, and add the newly detected aromatic substructures to the
        AromaticManager.

        :param molfile: the path to the molfile.
        :return: None.
        """
        cpd = compound.Compound.molfile_name(molfile)
        if cpd:
            aromatic_bonds = self.indigo_aromatic_bonds(molfile)
            cycles = cpd.separate_connected_components(aromatic_bonds)
            aromatic_substructures = self.construct_aromatic_entity(cpd, cycles)
            self.add_aromatic_substructures(aromatic_substructures)

    def indigo_aromatic_bonds(self, molfile: str) -> set:
        """
        Detect the aromatic bonds in the compound via Indigo method.

        :param molfile: the path to the molfile.
        :return: the set of aromatic bonds represented by first_atom_number and second_atom_number in the bond.
        """
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
    def fuse_cycles(cycles: list) -> list:
        """
        To fuse the cycles with shared atoms.

        :param cycles: the list of cycles represented by atom numbers.
        :return: the list of cleaned cycles.
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

    # def detect_aromatic_substructures_timeout(self, cpd: compound.Compound) -> None:
    #     """
    #     To detect the aromatic substructures in the compound and stop the search on timeout.
    #     :param cpd: the :class:`~md_harmonize.compound.Compound` entity.
    #     :return: None.
    #     """
    #
    #     try:
    #         with tools.timeout(seconds=200):
    #             self.detect_aromatic_substructures(cpd)
    #     except:
    #         print("Aromatic substructures in compound {0} can hardly be detected.".format(cpd.name))
    #         pass
    #     return

    def detect_aromatic_substructures_timeout(self, cpd: compound.Compound) -> None:
        """
        Detect the aromatic substructures in the compound and stop the search on timeout.

        :param cpd: the :class:`~md_harmonize.compound.Compound` entity.
        :return: None.
        """
        try:
            tools.timeout(self.detect_aromatic_substructures, (cpd,), seconds=200)
        except:
            print("Aromatic substructures in compound {0} can hardly be detected.".format(cpd.name))
            pass
        return

    def detect_aromatic_substructures(self, cpd: compound.Compound) -> None:
        """
        Detect all the aromatic substructures in the cpd, and update the bond type of aromatic bonds.

        :param cpd: the :class:`~md_harmonize.compound.Compound` entity.
        :return: None.
        """
        cpd.update_color_tuple()
        cycles = []
        for aromatic in self.aromatic_substructures:
            if all(aromatic.composition[key] <= cpd.composition[key] for key in aromatic.composition):
                mapping_matrix = BASS.make_mapping_matrix(aromatic, cpd, True, True, False)
                if mapping_matrix is not None:
                    
                    for assignment in BASS.find_mappings(aromatic.structure_matrix(resonance=False),
                                                         aromatic.distance_matrix,
                                                         cpd.structure_matrix(resonance=False),
                                                         cpd.distance_matrix, mapping_matrix):
                        cycle = set()
                        for bond in aromatic.bonds:
                            if aromatic.atoms[bond.first_atom_number].in_cycle and \
                                    aromatic.atoms[bond.second_atom_number].in_cycle:
                                first_atom_number = cpd.heavy_atoms[assignment[bond.first_atom_number]].atom_number
                                second_atom_number = cpd.heavy_atoms[assignment[bond.second_atom_number]].atom_number
                                cycle.add(first_atom_number)
                                cycle.add(second_atom_number)
                        cycles.append(cycle)
        fused_cycles = self.fuse_cycles(cycles)
        cpd.update_aromatic_bond_type(fused_cycles)
        return 

    @staticmethod
    def construct_aromatic_entity(cpd: compound.Compound, aromatic_cycles: list) -> list:
        """
        Construct the aromatic substructure entity based on the aromatic atoms.
        Here, we also include outside atoms that are connected to aromatic rings with double bonds.

        :param cpd: the :class:`~md_harmonize.compound.Compound` entity.
        :param aromatic_cycles: the list of aromatic cycles represented by atom numbers in the compound.
        :return: the list of constructed aromatic substructures.
        """
        count = 0
        aromatic_substructures = []
        for cycle in aromatic_cycles:
            bonds = [bond.clone() for bond in cpd.bonds if bond.first_atom_number in cycle and bond.second_atom_number
                     in cycle]
            seen_atoms = set(cycle)
            for atom_index in cycle:
                atom = cpd.atoms[atom_index]
                for neighbor_index in atom.neighbors:
                    connecting_bond = cpd.bond_lookup[ (atom_index, neighbor_index)]
                    if neighbor_index not in cycle:
                        # add outside atoms connecting to aromatic ring with double bond.
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
            aromatic_substructures.append(compound.Compound(cpd.name + "_" + str(count), atoms, bonds))
            count += 1
        return aromatic_substructures

    def extract_aromatic_substructures(self, cpd: compound.Compound) -> list:
        """
        Detect the aromatic substructures in a compound based on the aromatic atoms. This only applies to KEGG kcf file.

        :param cpd: the :class:`~md_harmonize.compound.Compound` entity.
        :return: the list of aromatic cycles represented by atom numbers.
        """
        aromatic_elements = ["C", "N"]
        aromatic_types = ["C8x", "C8y", "N4x", "N4y", "N5x", "N5y"]
        aromatic_cycles = []
        for cycle in cpd.cycles:
            aromatic_atoms = [cpd.atoms[index] for index in cycle if cpd.atoms[index].default_symbol in
                              aromatic_elements]
            is_aromatic = [atom.kat in aromatic_types for atom in aromatic_atoms].count(True) == len(aromatic_atoms)
            oxygen_count = [cpd.atoms[index].default_symbol for index in cycle].count("O")
            # see the middle in KEGG compound C03861.
            if not is_aromatic or (oxygen_count >= 2 and len(cycle) >= 6):
                continue
            aromatic_cycles.append(cycle)
        return self.fuse_cycles(aromatic_cycles)

