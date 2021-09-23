#!/usr/bin/python3

import collections

def kegg_reaction_parser(reaction):

    """
    This is to parse KEGG reaction file to a dictionary.

    eg:
    ENTRY       R00259                      Reaction
    NAME        acetyl-CoA:L-glutamate N-acetyltransferase
    DEFINITION  Acetyl-CoA + L-Glutamate <=> CoA + N-Acetyl-L-glutamate
    EQUATION    C00024 + C00025 <=> C00010 + C00624
    RCLASS      RC00004  C00010_C00024
                RC00064  C00025_C00624
    ENZYME      2.3.1.1
    PATHWAY     rn00220  Arginine biosynthesis
                rn01100  Metabolic pathways
                rn01110  Biosynthesis of secondary metabolites
                rn01210  2-Oxocarboxylic acid metabolism
                rn01230  Biosynthesis of amino acids
    MODULE      M00028  Ornithine biosynthesis, glutamate => ornithine
                M00845  Arginine biosynthesis, glutamate => acetylcitrulline => arginine
    ORTHOLOGY   K00618  amino-acid N-acetyltransferase [EC:2.3.1.1]
                K00619  amino-acid N-acetyltransferase [EC:2.3.1.1]
                K00620  glutamate N-acetyltransferase / amino-acid N-acetyltransferase [EC:2.3.1.35 2.3.1.1]
                K11067  N-acetylglutamate synthase [EC:2.3.1.1]
                K14681  argininosuccinate lyase / amino-acid N-acetyltransferase [EC:4.3.2.1 2.3.1.1]
                K14682  amino-acid N-acetyltransferase [EC:2.3.1.1]
                K22476  N-acetylglutamate synthase [EC:2.3.1.1]
                K22477  N-acetylglutamate synthase [EC:2.3.1.1]
                K22478  bifunctional N-acetylglutamate synthase/kinase [EC:2.3.1.1 2.7.2.8]
    DBLINKS     RHEA: 24295
    ///

    :param str reaction: the KEGG reaction description.
    :return: the dictionary of the parsed KEGG reaction.
    """
    reaction_dict = collections.defaultdict(list)
    key = ""
    for line in reaction:
        if line.startswith("  "):
            reaction_dict[key].append(line)
        else:




def kegg_kcf_parser(kcf):

    """
    This is to parse KEGG kcf file to a dictionary.

    eg:
    ENTRY       C00013                      Compound
    ATOM        9
                1   P1b P    22.2269  -20.0662
                2   O2c O    23.5190  -20.0779
                3   O1c O    21.0165  -20.0779
                4   O1c O    22.2851  -21.4754
                5   O1c O    22.2617  -18.4642
                6   P1b P    24.8933  -20.0837
                7   O1c O    24.9401  -21.4811
                8   O1c O    26.1797  -20.0662
                9   O1c O    24.9107  -18.4582
    BOND        8
                1     1   2 1
                2     1   3 1
                3     1   4 1
                4     1   5 2
                5     2   6 1
                6     6   7 1
                7     6   8 1
                8     6   9 2
    ///

    :param kcf: the kcf text
    :return: the dictionary of parsed kcf file.
    """
    compound_name, atom_count, atoms, bond_count, bonds = "", 0, [], 0, []
    for line in kcf:
        tokens = [item for item in line.split(' ') if item != '']
        if tokens[0] == "ENTRY":
            compound_name = tokens[1]
        elif tokens[0] == "ATOM":
            atom_count = tokens[1]
        elif tokens[0] == "BOND":
            bond_count = tokens[1]
        elif tokens[0] == "///":
            continue
        else:
            if len(atoms) < atom_count:
                atoms.append({"atom_symbol": tokens[2], "atom_number": int(tokens[0])-1, "x": tokens[3], "y": tokens[4], "kat": tokens[1]})
            else:
                bonds.append({"first_atom_number":tokens[1], "second_atom_number": tokens[2], "bond_type": tokens[3]})
    return {"compound_name": compound_name, "atoms": atoms, "bonds": bonds}

def kegg_rclass_parser(rclass):

    """
    This is to parse KEGG Rclass text file to a dictionary.

    eg:
    ENTRY       RC00001                     RClass
    DEFINITION  C1x-C8x:*-*:C2x+C2y-C8x+C8y
                N1y-N5y:*-*:C1y+C2x+C2x-C1y+C8x+C8x
    RPAIR       C00003_C00004    C00005_C00006
    REACTION    R00090 R00112 R00267 R00282 R01528 R02034 R02163 R02736
                R07141 R07171 R07172 R10517 R10518 R12579
    ENZYME      1.1.1.42        1.1.1.43        1.1.1.44        1.1.1.49
                1.1.1.285       1.1.1.351       1.1.1.363       1.6.1.1
                1.6.1.2         1.6.1.3         1.6.3.1         1.6.3.2
                1.6.3.3         1.6.3.4         1.6.5.9         1.6.99.1
                1.11.1.1        1.11.1.26
    PATHWAY     rn00480  Glutathione metabolism
    ORTHOLOGY   K00031  isocitrate dehydrogenase [EC:1.1.1.42]
                K00032  phosphogluconate 2-dehydrogenase [EC:1.1.1.43]
                K00033  6-phosphogluconate dehydrogenase [EC:1.1.1.44 1.1.1.343]
                K00036  glucose-6-phosphate 1-dehydrogenase [EC:1.1.1.49 1.1.1.363]
                K00322  NAD(P) transhydrogenase [EC:1.6.1.1]
                K00323  H+-translocating NAD(P) transhydrogenase [EC:1.6.1.2 7.1.1.1]
                K00324  H+-translocating NAD(P) transhydrogenase subunit alpha [EC:1.6.1.2 7.1.1.1]
                K00325  H+-translocating NAD(P) transhydrogenase subunit beta [EC:1.6.1.2 7.1.1.1]
                K00354  NADPH2 dehydrogenase [EC:1.6.99.1]
                K13937  hexose-6-phosphate dehydrogenase [EC:1.1.1.47 3.1.1.31]
    ///

    :param rclass: the rclass text description.
    :return: the dictionary of the parsed rclass.
    """



def rpair_parser(rclass, one_compound, the_other_compound):

    """
    This is to get one to one atom mappings between two compounds based on the rclass definition.
    :param list rclass: the list of rclass definition.
    :param one_compound:
    :type one_compound:
    :param the_other_compound:
    :type
    :return:
    """


