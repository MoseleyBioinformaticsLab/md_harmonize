#!/usr/bin/python3

import tools
import compound
import BASS


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

def extract_aromatic_substructures(cpd):
    """
    To detect the aromatic substructures iin a compound based on the aromatic types.
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
    return fuse_cycles(aromatic_cycles)

def construct_aromatic_entities(cpd, aromatic_cycles):
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

def remove_duplicate_aromatic_substructure(aromatic_substructures):
    """
    To remove duplicate aromatic substructures.
    :param aromatic_substructures:
    :return:
    """
    skip = set()
    unique = []
    n = len(aromatic_substructures)
    for i in range(n):
        if i not in skip:
            aromatic_1 = aromatic_substructures[i]
            for j in range(i+1, n):
                if j not in skip:
                    aromatic_2 = aromatic_substructures[j]
                    if len(aromatic_1.atoms) == len(aromatic_2.atoms) and len(aromatic_1.bonds) == len(aromatic_2.bonds):
                        mapping_matrix = BASS.
                        if mapping_matrix:
                            isomorphs = BASS.
                            if isomorphs != [] and isomorphs is not None:
                                skip.add(j)
            unique.append(aromatic_1)
    return unique

def detect_aromatic_substructures(cpd, aromatic_substructures):
    """
    Detect the aromatic substructures of the cpd.
    :param cpd:
    :param aromatic_substructures:
    :return:
    """

    aromatic_cycles = []
    for aromatic in aromatic_substructures:
        if aromatic.composition.issubset(cpd.composition):
            mmat = BASS.mappping_matrix(aromatic, cpd, True, False, False, True)
            if mmat:
                for
    return fuse_cycles(aromatic_cycles)

def construct_aromatic_substructure_set(kegg_compounds, save_file):
    """
    To construct the aromatic substructure set based on the KEGG atom types.
    :param kegg_kcf_directory:
    :param save_file:
    :return:
    """
    aromatic_substructures = []
    for cpd in kegg_compounds:
        aromatic_cycles = extract_aromatic_substructures(cpd)
        aromatic_substructures.extend(construct_aromatic_entities(cpd, aromatic_cycles))
    unique_aromatic_substructures = remove_duplicate_aromatic_substructure(aromatic_substructures)
    tools.save_to_jsonpickle(unique_aromatic_substructures, save_file)
    return unique_aromatic_substructures






