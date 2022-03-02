#!/usr/bin/python3
"""
MDH command-line interface

Usage:
    MDH -h | --help
    MDH --version
    MDH download <database_names> <working_directory>
    MDH standardize <database_names> <working_directory>
    MDH aromatize <database_names> <working_directory> <save_file> [--aromatic_manager=<aromatic_manager_file>]
    MDH initialize <database_names> <working_directory> <aromatic_manager_file> [--kegg_atom_mappings=<kegg_atom_mappings_file>]
    MDH harmonize <database_names> <working_directory>
    MDH test
    MDH test1
    MDH test2
    MDH test3    
    MDH test4 <database_names> <working_directory>

Options:
    -h, --help          Show this screen.
    --version           Show version.

"""
import glob
import os
import multiprocessing
import time
import collections
from . import KEGG_database_scraper
from . import KEGG_parser
from . import MetaCyc_parser
from . import aromatics
from . import tools
from . import compound
from . import openbabel_utils
from . import harmonization

scraper_dict = {"KEGG": KEGG_database_scraper}
parser_dict = {"KEGG": KEGG_parser, "MetaCyc": MetaCyc_parser}

def construct_cpd(molfile):
    if os.path.getsize(molfile) != 0:
        return compound.Compound.create(molfile)
    return None

def construct_kcf(file):
    if os.path.getsize(file) != 0:
        return parser_dict['KEGG'].create_compound_kcf(file)
    return None

def compound_construct_multiprocess(entities, function):
    with multiprocessing.Pool() as pool:
        results = pool.map(function, entities)
    return {cpd.compound_name: cpd for cpd in results if cpd}

def construct_cpd_components(compound_components):
    return compound.Compound(compound_components[0], compound_components[1], compound_components[2])


def atom_order_check(compound_dict1, compound_dict2):
    # the compound name for kcf compound does not contain "cpd:".
    for compound_name in compound_dict1:
        compound1 = compound_dict1[compound_name]
        if compound_name not in compound_dict2 and "cpd:"+compound_name not in compound_dict2:
            continue
        compound2 = compound_dict2[compound_name] if compound_name in compound_dict2 else compound_dict2["cpd:"+compound_name]
        k = min(len(compound1.atoms), len(compound2.atoms))
        for i in range(k):
            if compound1.atoms[i].x != compound2.atoms[i].x or compound1.atoms[i].y != compound2.atoms[i].y:
                print("{0} file has changed!".format(compound_name))
                break
    return

def KEGG_atom_mapping_correction(atom_index_mappings, atom_mappings):
    
    new_atom_mappings = {}
    for name in atom_mappings:
        rlcass, cpd1, cpd2 = name.split("_")
        if cpd1 in atom_index_mappings and cpd2 in atom_index_mappings:
            this_mapping = []
            for entry in atom_mappings[name]:
                (from_cpd, from_idx), (to_cpd, to_idx) = entry
                this_mapping.append(((from_cpd, atom_index_mappings[from_cpd][from_idx]), (to_cpd, atom_index_mappings[to_cpd][to_idx])))
            new_atom_mappings[name] = this_mapping
    return new_atom_mappings

def KEGG_atom_index_mapping(kcf_compounds, mol_compounds):

    atom_index_mappings = {}
    for cpd_name in kcf_compounds:
        kcf_cpd = kcf_compounds[cpd_name]
        if "cpd:" + cpd_name in mol_compounds:
            mol_cpd = mol_compounds["cpd:" + cpd_name]
            mappings = {}
            i, j = 0, 0
            while i < len(kcf_cpd.atoms):
                while j < len(mol_cpd.atoms) and (kcf_cpd.atoms[i].x != mol_cpd.atoms[j].x or kcf_cpd.atoms[i].y != mol_cpd.atoms[j].y):
                    j += 1
                mappings[i] = j
                i += 1
                j += 1
            if i != len(kcf_cpd.atoms):
                print("KEGG atom index mappings for compound {0} are not complete!!!".format(cpd_name))
            else:
                atom_index_mappings[cpd_name] = mappings
    return atom_index_mappings

def cli(args):

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
        # here to search the database directory, find the molfile subdirectory and create standardized_molfile subdirectory
        # and store the standardized molfile in it.
        database_names = args['<database_names>'].split(",")
        working_directory = args['<working_directory>']
        for database_name in database_names:
            from_path = working_directory + "sources/{0}/molfile".format(database_name)
            if not os.path.exists(from_path):
                raise OSError("The directory {0} does not exist.".format(from_path))
            to_path = working_directory + "/standardized/{0}/molfile".format(database_name)
            os.makedirs(to_path, exist_ok=True)
            molfiles = glob.glob(from_path + "/*")
            for molfile in molfiles:
                if os.path.getsize(molfile) != 0:
                    openbabel_utils.standardize_molfile(molfile, to_path)

    elif args['aromatize']:
        database_names = args['<database_names>'].split(",")
        working_directory = args['<working_directory>']
        save_file = args['<save_file>']
        if args['--aromatic_manager']:
            aromatic_manager = aromatics.AromaticManager.decode(tools.open_jsonpickle(args['--aromatic_manager']))
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
        tools.save_to_jsonpickle(aromatic_manager.encode(), save_file)
    
    elif args['initialize']:
        database_names = args['<database_names>'].split(",")
        working_directory = args['<working_directory>']
        aromatic_manager = aromatics.AromaticManager.decode(tools.open_jsonpickle(args['<aromatic_manager_file>']))
        print("count of aromatic substructures,", len(aromatic_manager.aromatic_substructures))
        to_directory = working_directory + "/initialized"
        for database_name in database_names:
            from_path =  working_directory + "standardized/{0}/molfile".format(database_name)
            if not os.path.exists(from_path):
                raise OSError("The directory {0} does not exist.".format(from_path))
            parser = parser_dict.get(database_name, None)
            # construct compounds
            compound_dict = {}
            molfiles = glob.glob(from_path + "/*")
            compound_dict = compound_construct_multiprocess(molfiles, construct_cpd)
        
            reaction_list = []
            # construct reactions
            if database_name == "KEGG":
                reaction_directory = working_directory + "/sources/KEGG/reaction/"
                rclass_directory = working_directory + "/sources/KEGG/rclass/"
                if not os.path.exists(reaction_directory):
                    raise OSError("The directory {0} does not exist.".format(reaction_directory))
                if not os.path.exists(rclass_directory):
                    raise OSError("The directory {0} does not exist.".format(rclass_directory))
                # construct KEGG compound via KEGG kcf files. This is used for atom mappings 
                kcf_files = glob.glob(working_directory + "/sources/KEGG/kcf/*")
                kcf_compounds = compound_construct_multiprocess(kcf_files, construct_kcf)
                
                original_files = glob.glob(working_directory + "/sources/KEGG/molfile/*")
                original_compounds = compound_construct_multiprocess(original_files, construct_cpd)
            
                # sanity check if the atom orders in the kcf file, original molfile and standardized molfile are the same.
                # two comparison: kcf - original; original - standardized.
                # for kcf - original: some molfiles have H atoms but the corresponding kcf files don't, and these Hs are in the middle, which can lead to the mess of atom mappings.
                atom_order_check(kcf_compounds, original_compounds)
                atom_order_check(original_compounds, compound_dict)
                if args['--kegg_atom_mappings']:
                    atom_mappings = tools.open_json(args['--kegg_atom_mappings'])
                else:
                    atom_mappings = parser.create_atom_mappings(rclass_directory, kcf_compounds)
                
                # from kcf atom mappings to mol file atom mappings.
                atom_mappings = KEGG_atom_mapping_correction(KEGG_atom_index_mapping(kcf_compounds, compound_dict), atom_mappings)
                reaction_list = parser.create_reactions(reaction_directory, compound_dict, atom_mappings)

            elif database_name == "MetaCyc":
                reaction_file = working_directory + "/sources/MetaCyc/reactions.dat"
                atom_mapping_file = working_directory + "/sources/MetaCyc/atom-mapping.dat"
                if not os.path.exists(reaction_file):
                    raise OSError("The file {0} does not exist.".format(reaction_file))
                if not os.path.exists(atom_mapping_file):
                    raise OSError("The file {0} does not exist.".format(atom_mapping_file))
                reaction_list = parser.create_reactions(reaction_file, atom_mapping_file, compound_dict)

                # update compound
            for cpd_name in compound_dict:
                cpd = compound_dict[cpd_name]
                # aromatic substructure dection
                print("aromatic substructure dection: ", cpd_name)
                aromatic_manager.detect_aromatic_substructures(cpd)
                cpd.define_bond_stereochemistry()
                cpd.curate_invalid_n()

            save_directory = to_directory + "/{0}".format(database_name)
            os.makedirs(save_directory, exist_ok=True)
            tools.save_to_jsonpickle({name: compound_dict[name].encode() for name in compound_dict}, save_directory + "/compounds.json")
            tools.save_to_jsonpickle(reaction_list, save_directory + "/reactions.json")

    elif args['harmonize']:
        # create harmonized compound manager first, then reaction harmonization manager.

        database_names = args['<database_names>'].split(",")
        working_directory = args['<working_directory>']
        compound_dict = []
        reaction_list = []
        for database_name in database_names:
            from_directory = working_directory + "/initialized/{0}/".format(database_name)
            if not os.path.exists(from_directory + "compounds.json"):
                raise OSError("The file {0} does not exist.".format(from_directory + "compounds.json"))
            if not os.path.exists(from_directory + "reactions.json"):
                raise OSError("The file {0} does not exist.".format(from_directory + "reactions.json"))
            compounds = tools.open_jsonpickle(from_directory + "compounds.json")
            compound_dict.append({name: compound.Compound(compounds[name][0], compounds[name][1], compounds[name][2]) for name in compounds})
            reaction_list.append(tools.open_jsonpickle(from_directory + "reactions.json"))
        print("start compound harmonization")
        compound_harmonization_manager = harmonization.harmonize_compound_list(compound_dict)
        # while we do reaction harmonization, we need to pay attention to compound harmonization without same structural representations.
        # this includes: R group, linear-circular-transformation, resonance.
        print("start reaction harmonization")
        reaction_harmonization_manager = harmonization.harmonize_reaction_list(reaction_list, compound_harmonization_manager)
        save_directory = working_directory + "/harmonized"
        os.makedirs(save_directory, exist_ok=True)
        tools.save_to_jsonpickle(reaction_harmonization_manager, save_directory + "/{0}.json".format("_".join(database_names)))

    elif args["test2"]:
        kcf_files = glob.glob("/mlab/data/hji236/projects/MDH_test/sources/KEGG/kcf/*")
        kcf_compounds = {}
        parser = parser_dict["KEGG"]
        for kcf_file in kcf_files:
            if os.path.getsize(kcf_file) != 0:
                this_compound = parser.create_compound_kcf(kcf_file)
                this_compound.color_compound()
                kcf_compounds[this_compound.compound_name] = this_compound
        filename_new = "/mlab/data/hji236/projects/MDH_test/kegg_atom_mappings_test.json"
        filename_pre = "/mlab/data/hji236/projects/reaction_harmonization/figure_repo/results/atom_mappings_KEGG_RCLASS.json"
        new_data = tools.open_json(filename_new)
        pre_data = tools.open_json(filename_pre)
        new_dict = collections.defaultdict(dict)
        pre_dict = collections.defaultdict(dict)
        for name in new_data:
            for item in new_data[name]:
                (from_cpd, from_idx), (to_cpd, to_idx) = item
                new_dict[name]["{0}-{1}".format(from_cpd, from_idx)] = "{0}-{1}".format(to_cpd, to_idx)    
        for name in pre_data:
            rclass, cpd1, cpd2 = name.split("_")
            name1 = "_".join([rclass, cpd2, cpd1])
            for idx1 in pre_data[name]:
                idx2 = pre_data[name][idx1]
                pre_dict[name1]["{0}-{1}".format(cpd1, idx1)] = "{0}-{1}".format(cpd2, idx2)
                pre_dict[name1]["{0}-{1}".format(cpd2, idx2)] = "{0}-{1}".format(cpd1, idx1)
                pre_dict[name]["{0}-{1}".format(cpd1, idx1)] = "{0}-{1}".format(cpd2, idx2)
                pre_dict[name]["{0}-{1}".format(cpd2, idx2)] = "{0}-{1}".format(cpd1, idx1)
        for name in new_dict:

            if name not in pre_dict:
                print(name, "not in the parsed previously")
                continue
            for key in new_dict[name]:
                value = new_dict[name][key]
                if key not in pre_dict[name] or pre_dict[name][key] != value:
                    cpd1, idx1 = key.split("-")
                    cpd2, idx2 = value.split("-")
                    cpd1 = kcf_compounds[cpd1]
                    cpd2 = kcf_compounds[cpd2]
                    color1 = cpd1.atoms[int(idx1)].color
                    color2 = cpd2.atoms[int(idx2)].color
                    if len(cpd1.color_groups()[color1]) == 1 and len(cpd2.color_groups()[color2]) == 1:
                        print(name, "is not parsed consistently")
                        break
                
    elif args["test4"]:
        database_names = args['<database_names>'].split(",")
        print(database_names)
        working_directory = args['<working_directory>']
        compound_dict = []                    
        for database_name in database_names:
            from_directory = working_directory + "/initialized/{0}/".format(database_name)
            if not os.path.exists(from_directory + "compounds.json"):
                raise OSError("The file {0} does not exist.".format(from_directory + "compounds.json"))
            compounds = tools.open_jsonpickle(from_directory + "compounds.json")
            print("construct compounds")
            compounds = compound_construct_multiprocess(list(compounds.values()), construct_cpd_components)
            compound_dict.append(compounds)
            print("complete compound construction")
        compound_harmonization_manager = harmonization.harmonize_compound_list(compound_dict)
        save_directory = working_directory + "/harmonized/"
        os.makedirs(save_directory, exist_ok=True)
        tools.save_to_jsonpickle(compound_harmonization_manager, working_directory + "{0}.json".format("_".join(database_names)))
    
    elif args["test3"]:
        compound_file = "/mlab/data/hji236/projects/MDH_test/initialized/MetaCyc/compounds.json"
        compounds = tools.open_json(compound_file)
        working_directory = "/mlab/data/hji236/projects/MDH_test/"
        reaction_file = working_directory + "/sources/MetaCyc/reactions.dat"
        atom_mapping_file = working_directory + "/sources/MetaCyc/atom-mapping.dat"
        parser = parser_dict["MetaCyc"]
        reaction_list = parser.create_reactions(reaction_file, atom_mapping_file, compounds)
        tools.save_to_jsonpickle(reaction_list, "/mlab/data/hji236/projects/MDH_test/initialized/MetaCyc/reactions.json")
        
    elif args["test1"]:
        kcf_files = glob.glob("/mlab/data/hji236/projects/MDH_test/sources/KEGG/kcf_test1/*")
        kcf_compounds = {}
        parser = parser_dict["KEGG"]
        for kcf_file in kcf_files:
            if os.path.getsize(kcf_file) != 0:
                this_compound = parser.create_compound_kcf(kcf_file)
                this_compound.color_compound()
                kcf_compounds[this_compound.compound_name] = this_compound
        for cpd in kcf_compounds:
            print(cpd, "atom number", len(kcf_compounds[cpd].atoms))
        rclass_directory = "/mlab/data/hji236/projects/MDH_test/sources/KEGG/rclass_test1/"
        atom_mappings = parser.create_atom_mappings(rclass_directory, kcf_compounds)
        to_file = "/mlab/data/hji236/projects/MDH_test/kegg_atom_mappings_test1.json"
        tools.save_to_json(atom_mappings, to_file)
        

    elif args['test']:
        aromatic_manager_file = "/mlab/data/hji236/projects/MDH_test/aromatic_manager.json"
        aromatic_manager = aromatics.AromaticManager.decode(tools.open_jsonpickle(aromatic_manager_file))
        kegg_molfiles = glob.glob("/mlab/data/hji236/projects/MDH_test/standardized/KEGG/molfile/*")
        kegg_dict = {}
        metacyc_molfiles = glob.glob("/mlab/data/hji236/projects/MDH_test/standardized/MetaCyc/molfile/*")
        metacyc_dict = {}

        print("construct kegg compounds")
        for molfile in kegg_molfiles:
            if os.path.getsize(molfile) != 0:
                try:
                    cpd = compound.Compound.create(molfile)
                    kegg_dict[cpd.compound_name] = cpd
                except Exception as e:
                    print(e)
        print("construct metacyc compounds")
        for molfile in metacyc_molfiles:
            if os.path.getsize(molfile) != 0:
                try:
                    cpd = compound.Compound.create(molfile)
                    metacyc_dict[cpd.compound_name] = cpd
                except Exception as e:
                    print(e)
        print("detect aromatic structure in kegg")
        for cpd_name in kegg_dict:
                cpd = kegg_dict[cpd_name]
                # aromatic substructure dection
                print("aromatic substructure dection", cpd.compound_name)
                aromatic_manager.detect_aromatic_substructures(cpd)
                cpd.define_bond_stereochemistry()
                cpd.curate_invalid_n()

        print("detect aromatic structure in metacyc")
        for cpd_name in metacyc_dict:
                cpd = metacyc_dict[cpd_name]
                # aromatic substructure dection
                print("aromatic substructure dection", cpd.compound_name)
                aromatic_manager.detect_aromatic_substructures(cpd)
                cpd.define_bond_stereochemistry()
                cpd.curate_invalid_n()

