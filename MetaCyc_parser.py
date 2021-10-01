#!/usr/bin/python3

"""
MetaCyc_parser.py is used to parse MetaCyc text data. 
Note: All MetaCyc reactions atom_mappings are stored in a single text file.
"""

import collections
import copy
import glob
from pathlib import Path
import compound
import ctfile
import tools
import reaction

def reaction_side_parser(reaction_side):
    """
    This is to parse FROM_SIDE or TO_SIDE in the reaction.

    eg: FROM-SIDE - (CPD-9147 0 8) (OXYGEN-MOLECULE 9 10)
	
    Information includes compound name and the start and end mapping atoms in this compound.
    The order of the atoms are the orders in the compound molfile.
    """
    i = 0
    compounds = collections.defaultdict(list)
    while i < len(reaction_side):
        if reaction_side[i] == "(":
            i += 1
            count = 1
            start_point = i
            while count > 0:
                if reaction_side[i] == "(":
                    count += 1
                if reaction_side[i] == ")":
                    count -= 1
                i += 1
            i -= 1
            inner_substring = reaction_side[start_point: i]
            if "(" in inner_substring:
                compound, n, s, e = inner_substring.split()
                compounds[compound[1:]].append((int(s), int(e)))
            else:
                compound, s, e = inner_substring.split()
                compounds[compound].append((int(s), int(e)))
        i += 1
    return compounds

def generate_one_to_one_mappings(from_side, to_side, indices):
    """

    :param from_side:
    :param to_side:
    :param indices:
    :return:
    """
    n = len(indices)
    from_index_dict = {}
    for cpd_name in from_side:
        start_index, end_index = from_side[cpd_name][0]
        for i in range(start_index, end_index+1):
            atom_index = i - start_index
            from_index_dict[i] = (cpd_name, atom_index)
    to_index_dict = {}
    for cpd_name in to_side:
        start_index, end_index = to_side[cpd_name][0]
        for i in range(start_index, end_index+1):
            atom_index = i - start_index
            to_index_dict[i] = (cpd_name, atom_index)
    one_to_one_mappings = []
    for to_index, from_index in enumerate(indices):
        one_to_one_mappings.append((from_index_dict[from_index], to_index_dict[to_index]))
    return one_to_one_mappings

def reaction_with_reaction_side_parser(atom_mappings):
    """
    This is to parse the MetaCyc reaction with atom mappings.

    eg:
    REACTION - RXN-11981
    NTH-ATOM-MAPPING - 1
    MAPPING-TYPE - NO-HYDROGEN-ENCODING
    FROM-SIDE - (CPD-12950 0 23) (WATER 24 24) 
    TO-SIDE - (CPD-12949 0 24) 
    INDICES - 0 1 2 3 5 4 7 6 9 10 11 13 12 14 15 16 17 8 18 19 21 20 22 24 23

    note: the INDICES are atom mappings between two sides of the reaction. 
    TO-SIDE[i] is mapped to FROM-SIDE[idx] for i, idx in enumerate(INDICES).
    Pay attention to the direction!

    """
    reaction_dicts = {}
    current_reaction = {}
    for line in atom_mappings:
        if line.startswith("#"):
            continue
        elif line.startswith("//"):
            reaction_dicts["ONE_TO_ONE_MAPPINGS"] = generate_one_to_one_mappings(current_reaction["FROM-SIDE"],current_reaction["TO-SIDE"], current_reaction["INDICES"])
            reaction_dicts[current_reaction['REACTION']] = copy.deepcopy(current_reaction)
            current_reaction = {}
        else:
            key, value = line.split(" - ")
            if key == "FROM-SIDE" or key == "TO-SIDE":
                current_reaction[key] = reaction_side_parser(value)
            else:
                current_reaction[key] = value
    return reaction_dicts

def reaction_parser(reaction_text):
    """
    This is used to parse MetaCyc reaction.

    eg:
    UNIQUE-ID - RXN-13583
    TYPES - Redox-Half-Reactions
    ATOM-MAPPINGS - (:NO-HYDROGEN-ENCODING (1 0 2) (((WATER 0 0) (HYDROXYLAMINE 1 2)) ((NITRITE 0 2))))
    CREDITS - SRI
    CREDITS - caspi
    IN-PATHWAY - HAONITRO-RXN
    LEFT - NITRITE
    ^COMPARTMENT - CCO-IN
    LEFT - PROTON
    ^COEFFICIENT - 5
    ^COMPARTMENT - CCO-IN
    LEFT - E-
    ^COEFFICIENT - 4
    ORPHAN? - :NO
    PHYSIOLOGICALLY-RELEVANT? - T
    REACTION-BALANCE-STATUS - :BALANCED
    REACTION-DIRECTION - LEFT-TO-RIGHT
    RIGHT - HYDROXYLAMINE
    ^COMPARTMENT - CCO-IN
    RIGHT - WATER
    ^COMPARTMENT - CCO-IN
    STD-REDUCTION-POTENTIAL - 0.1    
    //
    """
    reaction_dicts = {}
    current_reaction = collections.defaultdict(list)
    count_left, count_right, previous_key = 0, 0, ""
    for line in reaction_text:
        if line.startswith("#"):
            continue
        elif line.startswith("//"):
            while len(current_reaction["^COEFFICIENT"]) < count_left + count_right:
                current_reaction["^COEFFICIENT"].append(" ")
            while len(current_reaction["^COMPARTMENT"]) < count_left + count_right:
                current_reaction["^COMPARTMENT"].append(" ")
            reaction_dicts[current_reaction['UNIQUE-ID'][0]] = copy.deepcopy(current_reaction)
            current_reaction = collections.defaultdict(list)
            count_left, count_right, previous_key = 0, 0, ""
        elif line.startswith('/'):
            current_reaction[previous_key].append(line)
        else:
            key, value = line.split(" - ")
            if key == 'LEFT':
                count_left += 1
            if key == 'RIGHT':
                count_right += 1
            if key == "^COEFFICIENT" or key == "^COMPARTMENT":
                while len(current_reaction[key]) < count_left + count_right - 1:
                    current_reaction[key].append(" ")
                current_reaction[key].append(value)
                previous_key = key
    return reaction_dicts

def create_metacyc_compounds(compound_directory):
    """

    :param compound_directory:
    :return:
    """
    compound_files = glob.glob(compound_directory+"*")
    compounds = {}
    for compound_file in compound_files:
        compound_name = Path(compound_file).stem
        with open(compound_file, 'r') as infile:
            ct_object = ctfile.load(infile)
            compounds[compound_name] = compound.Compound.create(ct_object, compound_name)
    return compounds

def create_metacyc_reactions(reaction_file, atom_mapping_file, compounds):
    """
    To create MetaCyc reaction entities.
    :param reaction_file:
    :param atom_mapping_file:
    :param compounds:
    :return:
    """
    reaction_dict = reaction_parser(tools.open_text(reaction_file).split("\n"))
    atom_mappings = reaction_with_reaction_side_parser(tools.open_text(atom_mapping_file).split("\n"))
    reactions = []
    for reaction_name in reaction_dict:
        this_reaction = reaction_dict[reaction_name]
        coefficient_list = collections.deque(this_reaction["^COEFFICIENT"])
        coefficients = {}
        one_side_compounds = []
        if "LEFT" in this_reaction:
            for cpd_name in this_reaction["LEFT"]:
                cpd_coefficient = coefficient_list.popleft()
                one_side_compounds.append(compounds[cpd_name])
                coefficients[cpd_name] = int(cpd_coefficient) if cpd_coefficient != " " else 1
        the_other_side_compounds = []
        if "RIGHT" in this_reaction:
            for cpd_name in this_reaction["RIGHT"]:
                cpd_coefficient = coefficient_list.popleft()
                the_other_side_compounds.append(compounds[cpd_name])
                coefficients[cpd_name] = int(cpd_coefficient) if cpd_coefficient != " " else 1
        ecs = collections.defaultdict(list)
        if "EC-NUMBER" in this_reaction:
            for ec in this_reaction["EC-NUMBER"]:
                if "|" not in ec:
                    numbers = ec[3:].split(".")
                    ecs[len(numbers)].append(ec[3:])
                else:
                    numbers = ec[4:-1].split(".")
                    ecs[len(numbers)].append(ec[4:-1])
        reactions.append(reaction.Reaction(reaction_name, one_side_compounds, the_other_side_compounds, ecs,
                                           atom_mappings[reaction_name]["ONE_TO_ONE_MAPPINGS"], coefficients))
    return reactions












