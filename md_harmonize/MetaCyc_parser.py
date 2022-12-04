#!/usr/bin/python3

"""

md_harmonize.MetaCyc_parser
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This module provides functions to parse MetaCyc text data.

Note: All MetaCyc reactions atom_mappings are stored in a single text file.

"""

import collections
import copy

from . import tools
from . import reaction


def reaction_side_parser(reaction_side: str) -> dict:
    """
    This is to parse FROM_SIDE or TO_SIDE in the reaction.

    eg: FROM-SIDE - (CPD-9147 0 8) (OXYGEN-MOLECULE 9 10)

    Information includes compound name and the start and end atom index in this compound used for atom mappings.
    The order of the atoms are the orders in the compound molfile.

    :param reaction_side: the text description of reaction side.
    :return: the dictionary of compounds and the corresponding start and end atom index in the atom mappings.
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


def generate_one_to_one_mappings(from_side: dict, to_side: dict, indices: str) -> list:
    """
    To generate the one to one atom mappings between two the sides of a metabolic reaction.

    :param from_side: the dictionary of compounds with their corresponding start and end atom indices in the from_side.
    :param to_side: the dictionary of compounds with their corresponding start and end atom indices in the to_side.
    :param indices: the string representation of mapped atoms.
    :return: the list of mapped atoms between the two sides (from_index, to_index).
    """
    from_index_dict = {}
    mappings = [int(num) for num in indices.split(" ") if num != ""]
    
    for cpd_name in from_side:
        for (start_index, end_index) in from_side[cpd_name]:
            for i in range(start_index, end_index+1):
                atom_index = i - start_index
                from_index_dict[i] = (cpd_name, atom_index)
    to_index_dict = {}
    
    for cpd_name in to_side:
        for (start_index, end_index) in to_side[cpd_name]:
            for i in range(start_index, end_index+1):
                atom_index = i - start_index
                to_index_dict[i] = (cpd_name, atom_index)
    one_to_one_mappings = []
    
    for to_index, from_index in enumerate(mappings):
        one_to_one_mappings.append((from_index_dict[from_index], to_index_dict[to_index]))
    return one_to_one_mappings


def atom_mappings_parser(atom_mapping_text: list) -> dict:
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

    :param atom_mapping_text: the text descriptions of reactions with atom mappings.
    :return: the dictionary of reactions with atom mappings.
    """
    reaction_dicts = {}
    current_reaction = {}
    for line in atom_mapping_text:
        if line.startswith("#"):
            continue
        elif line.startswith("//"):
            reaction_dicts[current_reaction['REACTION']] = copy.deepcopy(current_reaction)
            reaction_dicts[current_reaction['REACTION']]["ONE_TO_ONE_MAPPINGS"] = \
                generate_one_to_one_mappings(current_reaction["FROM-SIDE"], current_reaction["TO-SIDE"],
                                             current_reaction["INDICES"])
            current_reaction = {}
        else:
            key = line.split(" - ")[0]
            value = " - ".join(line.split(" - ")[1:])
            if key == "FROM-SIDE" or key == "TO-SIDE":
                current_reaction[key] = reaction_side_parser(value)
            else:
                current_reaction[key] = value
    return reaction_dicts


def reaction_parser(reaction_text: list) -> dict:
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

    :param reaction_text: the text descriptions of MetaCyc reactions.
    :return: the dict of parsed MetaCyc reactions.
    """
    reaction_dicts = {}
    current_reaction = collections.defaultdict(list)
    count_left, count_right, previous_key = 0, 0, ""
    for line in reaction_text:
        if line.startswith("#"):
            continue
        if line.startswith("//"):
            while len(current_reaction["^COEFFICIENT"]) < count_left + count_right:
                current_reaction["^COEFFICIENT"].append(" ")
            while len(current_reaction["^COMPARTMENT"]) < count_left + count_right:
                current_reaction["^COMPARTMENT"].append(" ")
            reaction_dicts[current_reaction['UNIQUE-ID'][0]] = copy.deepcopy(current_reaction)
            # print(reaction_dicts[current_reaction['UNIQUE-ID'][0]])
            current_reaction = collections.defaultdict(list)
            count_left, count_right, previous_key = 0, 0, ""
            continue
        # if line.startswith('/'):
        #     current_reaction[previous_key].append(line)
        #     continue
        key = line.split(" - ")[0]
        value = " - ".join(line.split(" - ")[1:])
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


def create_reactions(reaction_file: str, atom_mapping_file: str) -> list:
    """
    To molfile_name MetaCyc reaction entities.

    :param reaction_file: the path to the reaction file.
    :param atom_mapping_file: the path to the atom mapping file.
    :return: the list of constructed :class:`~md_harmonize.reaction.Reaction` entities.
    """
    reaction_dict = reaction_parser(tools.open_text(reaction_file, encoding='cp1252').split("\n"))
    atom_mappings = atom_mappings_parser(tools.open_text(atom_mapping_file).split("\n"))
    reactions = []
    for reaction_name in reaction_dict:
        this_reaction = reaction_dict[reaction_name]
        coefficient_list = collections.deque(this_reaction["^COEFFICIENT"])
        coefficients = {}
        if "LEFT" in this_reaction:
            for cpd_name in this_reaction["LEFT"]:
                cpd_coefficient = coefficient_list.popleft()
                coefficients[cpd_name] = cpd_coefficient if cpd_coefficient != " " else "1"
        if "RIGHT" in this_reaction:
            for cpd_name in this_reaction["RIGHT"]:
                cpd_coefficient = coefficient_list.popleft()
                coefficients[cpd_name] = cpd_coefficient if cpd_coefficient != " " else "1"
        ecs = collections.defaultdict(list)
        if "EC-NUMBER" in this_reaction:
            for ec in this_reaction["EC-NUMBER"]:
                if "|" not in ec:
                    numbers = ec[3:].split(".")
                    ecs[len(numbers)].append(ec[3:])
                else:
                    numbers = ec[4:-1].split(".")
                    ecs[len(numbers)].append(ec[4:-1])
        this_mappings = atom_mappings[reaction_name]["ONE_TO_ONE_MAPPINGS"] if reaction_name in atom_mappings and \
                                                                               "ONE_TO_ONE_MAPPINGS" in \
                                                                               atom_mappings[reaction_name] else []
        reactions.append(reaction.Reaction(reaction_name, this_reaction["LEFT"], this_reaction["RIGHT"], ecs,
                                           this_mappings, coefficients))
    return reactions

