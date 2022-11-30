#!/usr/bin/python3
"""
MDH command-line interface

Usage:
    mdh -h | --help
    mdh --version
    mdh download <database_names> <working_directory>
    mdh standardize <database_names> <working_directory>
    mdh aromatize <database_names> <working_directory> <save_file> [--aromatic_manager=<aromatic_manager_file>] [--pickle]
    mdh initialize_compound <database_names> <working_directory> <aromatic_manager_file> [--parse_kegg_atom] [--pickle]
    mdh initialize_reaction <database_names> <working_directory> [--pickle]
    mdh harmonize_compound <database_names> <working_directory> [--pickle]
    mdh harmonize_reaction <database_names> <working_directory> [--pickle]
    mdh test4 <database_names> <working_directory> [--pickle]
    mdh test7 <k>
    mdh test1
    mdh test2
    mdh test3
    mdh test4

Options:
    -h, --help          Show this screen.
    --version           Show version.

"""
import glob
import os
import multiprocessing
from . import KEGG_database_scraper
from . import KEGG_parser
from . import MetaCyc_parser
from . import aromatics
from . import tools
from . import compound
from . import openbabel_utils
from . import harmonization
from typing import *
import faulthandler
faulthandler.enable()


scraper_dict = {"KEGG": KEGG_database_scraper}
parser_dict = {"KEGG": KEGG_parser, "MetaCyc": MetaCyc_parser}
save_functions = {"pickle": tools.save_to_pickle, "jsonpickle": tools.save_to_jsonpickle}
open_functions = {"pickle": tools.open_pickle, "jsonpickle": tools.open_jsonpickle}


def construct_compound_via_molfile(molfile: str) -> Optional[compound.Compound]:
    """
    To construct compound based on the molfile representation.

    :param molfile: the molfile representation of the compound
    :return: the :class:`~mdh.compound.Compound` entity.
    """
    if os.path.getsize(molfile) != 0:
        return compound.Compound.create(molfile)
    return None


def construct_compound_via_kcf(file: str) -> Optional[compound.Compound]:
    """
    To construct compound based on the kcf representation.

    :param file: the kcf representation of the compound.
    :return: the :class:`~mdh.compound.Compound` entity.
    """
    if os.path.getsize(file) != 0:
        return parser_dict['KEGG'].create_compound_kcf(file)
    return None


def construct_compound_via_components(compound_components: list) -> Optional[compound.Compound]:
    """
    To construct compound based on the compound components.

    :param compound_components: the compound components.
    :return: the :class:`~mdh.compound.Compound` entity.
    """
    try:
        this_compound = compound.Compound(compound_components[0], compound_components[1], compound_components[2])
    except:
        this_compound = None
        pass
    return this_compound


def compound_construct_all(entities: list, function) -> dict:
    """
    To construct compounds one by one.

    :param entities: the list of entities used for constructing compounds.
    :param function: the function used to construct compounds.
    :return: the dict of compounds.
    """
    compounds = {}
    for entity in entities:
        # print(entity)
        this_compound = function(entity)
        if this_compound:
            compounds[this_compound.compound_name] = this_compound
    return compounds


def compound_construct_multiprocess(entities: list, function) -> dict:
    """
    To construct compounds with multiprocessing.

    :param entities: the list of entities used for constructing compounds.
    :param function: the function used to construct compounds.
    :return: the dict of compounds.
    """
    with multiprocessing.Pool() as pool:
        results = pool.map(function, entities)
    return {cpd.compound_name: cpd for cpd in results if cpd}


def atom_order_check(compound_dict1: dict, compound_dict2:dict) -> None:
    """
    Here, we want to check if the order of the heavy atoms in the KEGG molfile and kcf representations is the same.
    Quality check.

    :param compound_dict1: one dict of compounds
    :param compound_dict2: the other dict of compounds.
    :return: None
    """
    # the compound name for kcf compound does not contain "cpd:".
    count = 0
    for compound_name in compound_dict1:
        compound1 = compound_dict1[compound_name]
        if compound_name not in compound_dict2 and "cpd:"+compound_name not in compound_dict2:
            continue
        compound2 = compound_dict2[compound_name] if compound_name in compound_dict2 else \
            compound_dict2["cpd:"+compound_name]
        atoms1, atoms2 = compound1.heavy_atoms, compound2.heavy_atoms
        if len(atoms1) != len(atoms2):
            count += 1
            continue
        for atom1, atom2 in zip(atoms1, atoms2):
            if atom1.x != atom2.x or atom1.y != atom2.y:
                count += 1
                print("{0} file has changed!".format(compound_name))
                break
    print(count)
    return None


def KEGG_atom_mapping_correction(atom_index_mappings: dict, atom_mappings: dict) -> dict:
    """
    The corrected atom mappings between KEGG compounds composed of molfile representations.

    :param atom_index_mappings: the dict of atom index mappings between kcf and molfile representations.
    :param atom_mappings: the dict of atom mappings between compounds derived from kcf representations.
    :return: the dict of atom mappings corrected for molfile representations.
    """
    new_atom_mappings = {}
    for name in atom_mappings:
        rlcass, cpd1, cpd2 = name.split("_")
        if cpd1 in atom_index_mappings and cpd2 in atom_index_mappings:
            this_mapping = []
            for entry in atom_mappings[name]:
                (from_cpd, from_idx), (to_cpd, to_idx) = entry
                this_mapping.append(((from_cpd, atom_index_mappings[from_cpd][from_idx]),
                                     (to_cpd, atom_index_mappings[to_cpd][to_idx])))
            new_atom_mappings[name] = this_mapping
    return new_atom_mappings


def KEGG_atom_index_mapping(kcf_compounds: dict, mol_compounds: dict) -> dict:
    """
    To map the atom index between kcf and molfile representations for KEGG compounds.

    :param kcf_compounds: the dict of compounds composed of kcf representations.
    :param mol_compounds:the dict of compounds composed of molfile representations.
    :return: the dict of atom index mappings.
    """
    atom_index_mappings = {}
    for cpd_name in kcf_compounds:
        kcf_cpd = kcf_compounds[cpd_name]
        if "cpd:" + cpd_name in mol_compounds:
            mol_cpd = mol_compounds["cpd:" + cpd_name]
            mappings = {}
            i, j = 0, 0
            while i < len(kcf_cpd.atoms):
                while j < len(mol_cpd.atoms) and (kcf_cpd.atoms[i].x != mol_cpd.atoms[j].x or
                                                  kcf_cpd.atoms[i].y != mol_cpd.atoms[j].y):
                    j += 1
                mappings[i] = j
                i += 1
                j += 1
            if i != len(kcf_cpd.atoms):
                print("KEGG atom index mappings for compound {0} are not complete!!!".format(cpd_name))
            else:
                atom_index_mappings[cpd_name] = mappings
    return atom_index_mappings


def parse_reactions(compounds: dict, reactions: list) -> list:
    """
    To parse the reaction by linking to the compound entities.

    :param compounds: the dict of compound entities.
    :param reactions: the list of reaction entities.
    :return: the list of reaction entities.
    """

    for reaction in reactions:
        reaction.one_side = [compounds[name] for name in reaction.one_side if name in compounds]
        reaction.the_other_side = [compounds[name] for name in reaction.the_other_side if name in compounds]
    return reactions


def cli(args):

    save_function = save_functions["jsonpickle"]
    open_function = open_functions["jsonpickle"]

    if args['--pickle']:
        save_function = save_functions["pickle"]
        open_function = open_functions["pickle"]

    if args['download']:
        database_names = args['<database_names>'].split(",")
        working_directory = args['<working_directory>']
        for database_name in database_names:
            path = os.path.join(working_directory + "/sources", database_name)
            os.makedirs(path, exist_ok=True)
            try:
                scraper_dict[database_name].download(path)
            except KeyError:
                raise KeyError("Scraper for {0} database does not exits".format(database_name))

    elif args['standardize']:
        # here to search the database directory, find the molfile subdirectory and create standardized_molfile
        # subdirectory and store the standardized molfile in it.
        database_names = args['<database_names>'].split(",")
        working_directory = args['<working_directory>']
        for database_name in database_names:
            from_path = working_directory + "sources/{0}/molfile".format(database_name)
            if not os.path.exists(from_path):
                raise OSError("The directory {0} does not exist.".format(from_path))
            molfiles = glob.glob(from_path + "/*")
            if len(molfiles) > 25000:
                k = len(molfiles) // 25000
                for i in range(k+1):
                    to_path = working_directory + "/standardized/{0}/molfile_{1}".format(database_name, i)
                    os.makedirs(to_path, exist_ok=True)
                    for j in range(25000* i, min(25000* (i+1), len(molfiles))):
                        molfile = molfiles[j]
                        if os.path.getsize(molfile) != 0:
                            openbabel_utils.standardize_molfile(molfile, to_path)
            else:
                to_path = working_directory + "/standardized/{0}/molfile".format(database_name)
                os.makedirs(to_path, exist_ok=True)
                for molfile in molfiles:
                    if os.path.getsize(molfile) != 0:
                        openbabel_utils.standardize_molfile(molfile, to_path)

    elif args['aromatize']:
        database_names = args['<database_names>'].split(",")
        working_directory = args['<working_directory>']
        save_file = args['<save_file>']
        if args['--aromatic_manager']:
            aromatic_manager = aromatics.AromaticManager.decode(open_function(args['--aromatic_manager']))
        else:
            aromatic_manager = aromatics.AromaticManager()
        for database_name in database_names:
            from_path = working_directory + "standardized/{0}/molfile".format(database_name)
            if not os.path.exists(from_path):
                raise OSError("The directory {0} does not exist.".format(from_path))
            molfiles = glob.glob(from_path + "/*")
            for molfile in molfiles:
                if os.path.getsize(molfile) != 0:
                    aromatic_manager.indigo_aromatize(molfile)
        if "KEGG" in database_names:
            kcf_path = working_directory + "/sources/KEGG/kcf"
            if os.path.exists(kcf_path):
                kcf_files = glob.glob(kcf_path + "/*")
                parser = parser_dict['KEGG']
                for kcf_file in kcf_files:
                    if os.path.getsize(kcf_file) != 0:
                        kcf_cpd = parser.create_compound_kcf(kcf_file)
                        aromatic_manager.kegg_aromatize(kcf_cpd)
        save_function(aromatic_manager.encode(), save_file)
    
    elif args['initialize_compound']:
        database_names = args['<database_names>'].split(",")
        working_directory = args['<working_directory>']
        aromatic_manager = aromatics.AromaticManager.decode(open_function(args['<aromatic_manager_file>']))
        to_directory = working_directory + "/initialized"
        for database_name in database_names:
            if "_" in database_name:
                database_name, k = database_name.split("_")
                from_path = working_directory + "standardized/{0}/molfile_{1}".format(database_name, k)
                save_directory = to_directory + "/{0}".format(database_name)
                save_file = save_directory + "/compounds_{0}.json".format(k)
            else:
                from_path = working_directory + "standardized/{0}/molfile".format(database_name)
                save_directory = to_directory + "/{0}".format(database_name)
                save_file = save_directory + "/compounds.json"

            os.makedirs(save_directory, exist_ok=True)

            if not os.path.exists(from_path):
                raise OSError("The directory {0} does not exist.".format(from_path))
            parser = parser_dict.get(database_name, None)
            # construct compounds
            molfiles = glob.glob(from_path + "/*")
            compound_list = compound_construct_all(molfiles, construct_compound_via_molfile)

            if args['--parse_kegg_atom']:
                # I choose to parse the kegg atom mapings here, since the kegg atom mappings are only related to
                # compounds rather than reactions.
                # sanity check if the atom orders in the kcf file, original molfile and standardized molfile are the
                # same.
                # two comparison: kcf - original; original - standardized.
                # for kcf - original: some molfiles have H atoms but the corresponding kcf files don't, and these Hs are
                # in the middle, which can lead to the mess of atom mappings.

                rclass_directory = working_directory + "/sources/KEGG/rclass/"
                if not os.path.exists(rclass_directory):
                    raise OSError("The directory {0} does not exist.".format(rclass_directory))
                # construct KEGG compound via KEGG kcf files. This is used for atom mappings
                kcf_files = glob.glob(working_directory + "/sources/KEGG/kcf/*")
                kcf_compounds = compound_construct_multiprocess(kcf_files, construct_compound_via_kcf)
                original_files = glob.glob(working_directory + "/sources/KEGG/molfile/*")
                original_compounds = compound_construct_multiprocess(original_files, construct_compound_via_molfile)
                atom_order_check(kcf_compounds, original_compounds)
                atom_order_check(original_compounds, compound_list)
                atom_mappings = parser.create_atom_mappings(rclass_directory, kcf_compounds)
                # atom_mappings = tools.open_jsonpickle(working_directory + "/kegg_atom_mappings_test.json")
                save_function(atom_mappings, working_directory + "/kegg_atom_mappings.json")
                atom_mappings = KEGG_atom_mapping_correction(KEGG_atom_index_mapping(kcf_compounds, compound_list),
                                                             atom_mappings)
                save_function(atom_mappings, working_directory + "/kegg_atom_mappings_IC.json")

            for cpd_name in compound_list:
                cpd = compound_list[cpd_name]
                # aromatic substructure detection
                aromatic_manager.detect_aromatic_substructures(cpd)
                cpd.define_bond_stereochemistry()
                cpd.curate_invalid_n()

            save_function([compound_list[name].encode() for name in compound_list], save_file)

    elif args['initialize_reaction']:

        database_names = args['<database_names>'].split(",")
        working_directory = args['<working_directory>']
        to_directory = working_directory + "/initialized"
        for database_name in database_names:
            save_directory = to_directory + "/{0}".format(database_name)
            os.makedirs(save_directory, exist_ok=True)

            # if not os.path.exists(save_directory + "/compounds.json"):
            #     raise OSError("Please construct the {0} compounds first!".format(database_name))

            # compounds = open_function(save_directory + "/compounds.json")
            # compound_dict = compound_construct_all(compounds, construct_compound_via_components)

            reaction_list = []
            parser = parser_dict.get(database_name, None)

            if database_name == "KEGG":
                reaction_directory = working_directory + "/sources/KEGG/reaction/"
                if not os.path.exists(reaction_directory):
                    raise OSError("The directory {0} does not exist.".format(reaction_directory))
                if not os.path.exists(working_directory + "/kegg_atom_mappings_IC.json"):
                    raise OSError("The atom mappings of KEGG compounds have not been generated.")
                atom_mappings = open_function(working_directory + "/kegg_atom_mappings_IC.json")
                reaction_list = parser.create_reactions(reaction_directory, atom_mappings)

            elif database_name == "MetaCyc":
                reaction_file = working_directory + "/sources/MetaCyc/reactions.dat"
                atom_mapping_file = working_directory + "/sources/MetaCyc/atom-mapping.dat"
                if not os.path.exists(reaction_file):
                    raise OSError("The file {0} does not exist.".format(reaction_file))
                if not os.path.exists(atom_mapping_file):
                    raise OSError("The file {0} does not exist.".format(atom_mapping_file))
                reaction_list = parser.create_reactions(reaction_file, atom_mapping_file)
            save_function(reaction_list, save_directory + "/reactions.json")

    elif args['harmonize_compound']:
        # create harmonized compound manager first, then reaction harmonization manager.

        database_names = args['<database_names>'].split(",")
        working_directory = args['<working_directory>']
        compound_list = []
        # reaction_list = []
        save_directory = working_directory + "/harmonized"
        os.makedirs(save_directory, exist_ok=True)
        save_file = save_directory + "/{0}_initial_harmonized_compounds.json".format("_".join(database_names))

        if os.path.exists(save_file):
            print("we find the initial_harmonized_compounds")
            return

        for database_name in database_names:
            if "_" in database_name:
                database_name, k = database_name.split("_")
                from_directory = working_directory + "/initialized/{0}/".format(database_name)
                compound_file = from_directory + "compounds_{0}.json".format(k)
            else:
                from_directory = working_directory + "/initialized/{0}/".format(database_name)
                compound_file = from_directory + "compounds.json"
            if not os.path.exists(compound_file):
                raise OSError("The file {0} does not exist.".format(compound_file))

            compounds = open_function(compound_file)
            # this should be converted to compounds since it was saved in list format later.
            compound_parsed = compound_construct_multiprocess(compounds, construct_compound_via_components)
            compound_list.append(compound_parsed)
            # reactions = open_function(from_directory + "reactions.json")
            # reaction_parsed = parse_reactions(compound_parsed, reactions)
            # reaction_list.append(reaction_parsed)
        print("start compound harmonization")
        compound_harmonization_manager = harmonization.harmonize_compound_list(compound_list)
        initial_harmonized_compounds = compound_harmonization_manager.save_manager()
        save_function(initial_harmonized_compounds, save_file)

    elif args['harmonize_reaction']:

        database_names = args['<database_names>'].split(",")
        working_directory = args['<working_directory>']
        to_directory = working_directory + "/harmonized"

        initial_compound_pairs = to_directory + "/{0}_initial_harmonized_compounds.json".format("_".join(database_names))
        if not os.path.exists(initial_compound_pairs):
            raise OSError("Please do compound harmonization first.")

        compound_list = []
        reaction_list = []
        for database_name in database_names:

            compound_file = working_directory + "/initialized/{0}/compounds.json".format(database_name)
            if not os.path.exists(compound_file):
                raise OSError("Please construct {0} compounds first.".format(database_name))
            compounds = compound_construct_multiprocess(open_function(compound_file), construct_compound_via_components)
            compound_list.append(compounds)

            reaction_file = working_directory + "/initialized/{0}/reactions.json".format(database_name)
            if not os.path.exists(reaction_file):
                raise OSError("Please construct {0} reactions first.".format(database_name))
            reactions = open_function(reaction_file)
            parsed_reactions = parse_reactions(compounds, reactions)
            reaction_list.append(parsed_reactions)

        initial_harmonized_compounds = open_function(initial_compound_pairs)
        compound_harmonization_manager = harmonization.CompoundHarmonizationManager.create_manager(compound_list, initial_harmonized_compounds)

        # while we do reaction harmonization, we need to pay attention to compound harmonization without same structural
        # representations.
        # this includes: R group, linear-circular-transformation, resonance.

        print("start reaction harmonization")
        reaction_harmonization_manager = harmonization.harmonize_reaction_list(reaction_list, compound_harmonization_manager)
        # # let's just save the list of harmonized names first.
        harmonized_reactions = reaction_harmonization_manager.save_manager()
        save_function(harmonized_reactions, to_directory + "/{0}_harmonized_reactions.json".format("_".join(database_names)))

    elif args["test1"]:
        # kegg_compound_1_file = "/mlab/data/hji236/projects/MDH_test/sources/KEGG/kcf/cpd:C05670"
        # kegg_compound_2_file = "/mlab/data/hji236/projects/MDH_test/sources/KEGG/kcf/cpd:C06114"

        # kegg_compound_1_file = "/mlab/data/hji236/projects/MDH_test/sources/KEGG/kcf/cpd:C01063"
        # kegg_compound_2_file = "/mlab/data/hji236/projects/MDH_test/sources/KEGG/kcf/cpd:C09813"
        # kegg_1 = construct_compound_via_kcf(kegg_compound_1_file)
        # kegg_2 = construct_compound_via_kcf(kegg_compound_2_file)
        # compounds = {kegg_1.compound_name: kegg_1, kegg_2.compound_name: kegg_2}

        cpd_names = ["cpd:C00037", "cpd:C01921", "cpd:C00051", "cpd:C01921"]
        cpd = ["/mlab/data/hji236/projects/MDH_test/sources/KEGG/kcf/" + name for name in cpd_names]
        compounds = compound_construct_all(cpd, construct_compound_via_kcf)
        rclass_dir = "/mlab/data/hji236/projects/MDH_test/sources/KEGG/rclass_target/"
        atom_mappings = parser_dict["KEGG"].create_atom_mappings(rclass_dir, compounds)
        print(atom_mappings)

    elif args["test2"]:

        # test linear circular transformation.
        kegg_file = "/mlab/data/hji236/projects/MDH_test/standardized/KEGG/molfile/cpd:C00508.mol"
        metacyc_file = "/mlab/data/hji236/projects/MDH_test/standardized/MetaCyc/molfile/L-RIBULOSE.mol"

        kegg_cpd = construct_compound_via_molfile(kegg_file)
        print("right after construction", kegg_cpd.has_cycle)
        # kegg has cycle
        metacyc_cpd = construct_compound_via_molfile(metacyc_file)

        relationship1, atom_mappings1 = kegg_cpd.circular_pair_relationship(metacyc_cpd)
        # relationship2, atom_mappings2 = metacyc_cpd.circular_pair_relationship(kegg_cpd)

        print("relationship1, atom_mappings1", relationship1, atom_mappings1)
        # print("relationship2, atom_mappings2", relationship2, atom_mappings2)


    elif args["test3"]:

        kegg_file = "/mlab/data/hji236/projects/MDH_test/standardized/KEGG/molfile/cpd:C20416.mol"
        metacyc_file = "/mlab/data/hji236/projects/MDH_test/standardized/MetaCyc/molfile/CPD-12028.mol"
        kegg_cpd = construct_compound_via_molfile(kegg_file)
        metacyc_cpd = construct_compound_via_molfile(metacyc_file)
        kegg_cpd.color_compound(r_groups=True, bond_stereo=False, atom_stereo=False, resonance=False, isotope_resolved=False, charge=False)
        metacyc_cpd.color_compound(r_groups=True, bond_stereo=False, atom_stereo=False, resonance=False, isotope_resolved=False, charge=False)
        print("do resonant mappings")
        resonant_mappings = kegg_cpd.map_resonance(metacyc_cpd, r_distance=False)
        print(resonant_mappings)


    elif args["test4"]:

        kegg_miss_file = "/mlab/data/hji236/projects/MDH_results/KEGG_miss.json"
        metacyc_miss_file = "/mlab/data/hji236/projects/MDH_results/MetaCyc_miss.json"

        kegg_miss = tools.open_json(kegg_miss_file)
        metacyc_miss = tools.open_json(metacyc_miss_file)

        kegg_miss = {key: "cpd:" + value for (key, value) in kegg_miss}

        kegg_names = set(kegg_miss.values())
        metacyc_names = set(metacyc_miss.values())
        hmd_names = set(list(kegg_miss.keys()) + list(metacyc_miss.keys()))

        kegg_cpds = compound_construct_all(["/mlab/data/hji236/projects/MDH_test/standardized/KEGG/molfile/" + name + ".mol" for name in kegg_names], construct_compound_via_molfile)
        metacyc_cpds = compound_construct_all(["/mlab/data/hji236/projects/MDH_test/standardized/MetaCyc/molfile/" + name + ".mol" for name in metacyc_names], construct_compound_via_molfile)
        hmd_cpds = compound_construct_all(["/mlab/data/hji236/projects/MDH_test/standardized/HMD/molfile/" + name + ".mol" for name in hmd_names], construct_compound_via_molfile)

        kegg_no_structure = {}
        kegg_formula_issue = {}
        metacyc_no_structure = {}
        metacyc_formula_issue = {}

        for hmd in kegg_miss:
            if hmd not in hmd_cpds or kegg_miss[hmd] not in kegg_cpds:
                kegg_no_structure[hmd] = kegg_miss[hmd]
            elif hmd_cpds[hmd].formula() != kegg_cpds[kegg_miss[hmd]].formula():
                kegg_formula_issue[hmd] = kegg_miss[hmd]

        for hmd in metacyc_miss:
            if hmd not in hmd_cpds or metacyc_miss[hmd] not in metacyc_cpds:
                metacyc_no_structure[hmd] = metacyc_miss[hmd]
            elif hmd_cpds[hmd].formula() != metacyc_cpds[metacyc_miss[hmd]].formula():
                metacyc_formula_issue[hmd] = metacyc_miss[hmd]

        print("kegg no structure ", len(kegg_no_structure))
        print("kegg formula issue ", len(kegg_formula_issue))
        print("metacyc no structure ", len(metacyc_no_structure))
        print("metacyc formula issue ", len(metacyc_formula_issue))

        tools.save_to_json(kegg_no_structure, "/mlab/data/hji236/projects/MDH_results/KEGG_no_structure.json")
        tools.save_to_json(kegg_formula_issue, "/mlab/data/hji236/projects/MDH_results/KEGG_formula_issue.json")
        tools.save_to_json(metacyc_no_structure, "/mlab/data/hji236/projects/MDH_results/MetaCyc_no_structure.json")
        tools.save_to_json(metacyc_formula_issue, "/mlab/data/hji236/projects/MDH_results/MetaCyc_formula_issue.json")

    elif args["test5"]:
        database_name = args['<database_names>']
        working_directory = args['<working_directory>']
        # from_directory = working_directory + "/initialized/"

        # r group testing
        kegg_compound_file = "/scratch/hji236/MDH_test/standardized/KEGG/molfile/cpd:C01371.mol"
        metacyc_compound_file = "/scratch/hji236/MDH_test/standardized/MetaCyc/molfile/Alkanes.mol"


        #
        # kegg_compound_file = "/scratch/hji236/MDH_test/standardized/KEGG/molfile/cpd:C04618.mol"
        # metacyc_compound_file = "/scratch/hji236/MDH_test/standardized/MetaCyc/molfile/CPD-13230.mol"
        #
        kegg_compound = construct_compound_via_molfile(kegg_compound_file)
        metacyc_compound = construct_compound_via_molfile(metacyc_compound_file)
        print(kegg_compound.compound_name, metacyc_compound.compound_name)
        #
        # # kegg_compound = kegg_compound_parsed["cpd:C04618"]
        # # metacyc_compound = metacyc_compound_parsed["CPD-13230"]
        relationship, mapping = kegg_compound.with_r_pair_relationship(metacyc_compound)
        print(relationship, mapping)


        # resonance testing
        # kegg_compound_file = "/scratch/hji236/MDH_test/standardized/KEGG/molfile/cpd:C11821.mol"
        # metacyc_compound_file = "/scratch/hji236/MDH_test/standardized/MetaCyc/molfile/5-HYDROXYISOURATE.mol"
        # kegg_compound = construct_compound_via_molfile(kegg_compound_file)
        # metacyc_compound = construct_compound_via_molfile(metacyc_compound_file)
        # resonant_mappings = kegg_compound.map_resonance(metacyc_compound, r_distance=False)
        # if resonant_mappings:
        #     relationship, atom_mappings = kegg_compound.optimal_resonant_mapping(metacyc_compound, resonant_mappings)
        #     print(atom_mappings)
        #
        # # linear circular testing
        # kegg_compound_file = "/scratch/hji236/MDH_test/standardized/KEGG/molfile/cpd:C05345.mol"
        # metacyc_compound_file = "/scratch/hji236/MDH_test/standardized/MetaCyc/molfile/CPD-15709.mol"
        # kegg_compound = construct_compound_via_molfile(kegg_compound_file)
        # metacyc_compound = construct_compound_via_molfile(metacyc_compound_file)
        # relationship, atom_mappings = kegg_compound.circular_pair_relationship(metacyc_compound)
        # print(relationship, atom_mappings)



    # kegg_compound = construct_compound_via_molfile(kegg_compound_file)
    # metacyc_compound = construct_compound_via_molfile(metacyc_compound_file)
    # # compounds = tools.open_jsonpickle(from_directory + "KEGG/compounds.json")
    # # kegg_compound_parsed = { name: compound.Compound(compounds[name][0], compounds[name][1], compounds[name][2]) for
    # #                         name in compounds }
    # # compounds = tools.open_jsonpickle(from_directory + "MetaCyc/compounds.json")
    # # metacyc_compound_parsed = {name: compound.Compound(compounds[name][0], compounds[name][1], compounds[name][2]) for
    # #                         name in compounds}

    #     database_names = args['<database_names>'].split(",")
    #     ks = args['<ks>'].split(",")
    #     working_directory = args['<working_directory>']
    #     compound_dict = []
    #     for i, database_name in enumerate(database_names):
    #         from_directory = working_directory + "/initialized/{0}/".format(database_name)
    #         if ks[i] == "*":
    #             files = glob.glob("{0}compounds*.json".format(from_directory))
    #         else:
    #             files = glob.glob("{0}compounds_{1}.json".format(from_directory, ks[i]))
    #         if not files:
    #             raise OSError("The file {0} does not exist.".format(from_directory + "compounds.json"))
    #         compounds = {}
    #         for file in files:
    #             compounds.update(tools.open_jsonpickle(file))
    #         print("construct compounds")
    #         compounds = compound_construct_multiprocess(list(compounds.values()), construct_cpd_components)
    #         compound_dict.append(compounds)
    #         print("complete compound construction")
    #     compound_harmonization_manager = harmonization.harmonize_compound_list(compound_dict)
    #     save_directory = working_directory + "/harmonized/"
    #     os.makedirs(save_directory, exist_ok=True)
    #     tools.save_to_json(compound_harmonization_manager.save_manager(), working_directory +
    #                        "{0}_{1}.json".format("_".join(database_names), "_".join(ks)))
    #
    #
    # # construct compounds in MDH, sub group
    # elif args['test7']:
    #     aromatic_manager_file = "/scratch/hji236/MDH_test/aromatic_manager_1.json"
    #     aromatic_manager = aromatics.AromaticManager.decode(tools.open_jsonpickle(aromatic_manager_file))
    #     k = args['<k>']
    #     file_dir = "/scratch/hji236/MDH_test/standardized/HMD/molfile_split/sub_{0}".format(k)
    #     to_directory = "/scratch/hji236/MDH_test/initialized/HMD/"
    #     compound_dict = {}
    #     molfiles = glob.glob(file_dir + "/*")
    #     # compound_dict = compound_construct_multiprocess(molfiles, construct_cpd)
    #     for file in molfiles:
    #         print(file)
    #         compound = construct_cpd(file)
    #         compound_dict[compound.compound_name] = compound
    #     for cpd_name in compound_dict:
    #         cpd = compound_dict[cpd_name]
    #         # aromatic substructure detection
    #         print(cpd_name)
    #         aromatic_manager.detect_aromatic_substructures(cpd)
    #         cpd.define_bond_stereochemistry()
    #         cpd.curate_invalid_n()
    #     tools.save_to_jsonpickle({name: compound_dict[name].encode() for name in compound_dict},
    #                              to_directory + "/compounds_{0}.json".format(k))
    #

    
