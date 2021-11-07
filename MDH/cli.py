#!/usr/bin/python3
"""
MDH command-line interface

Usage:
    MDH -h | --help
    MDH --version
    MDH download <database_names> <working_directory>
    MDH standardize <database_names> <working_directory>
    MDH aromatize <database_names> <working_directory> <save_file> [--aromatic_manager=<aromatic_manager_file>]
    MDH initialize <database_names> <working_directory> <aromatic_manager_file>
    MDH harmonize <database_names> <working_directory>
    MDH test
    MDH test1

Options:
    -h, --help          Show this screen.
    --version           Show version.

"""
import glob
import os
import multiprocessing
import time
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

def function_multiprocess(entities, function):

    with multiprocessing.Pool() as pool:
        results = pool.map(function, entities)
    return {cpd.compound_name: cpd for cpd in results if cpd}


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
        to_directory = working_directory + "/initialized"
        for database_name in database_names:
            from_path =  working_directory + "standardized/{0}/molfile".format(database_name)
            if not os.path.exists(from_path):
                raise OSError("The directory {0} does not exist.".format(from_path))
            parser = parser_dict[database_name]
            # construct compounds
            compound_dict = {}
            molfiles = glob.glob(from_path + "/*")
            for molfile in molfiles:
                if os.path.getsize(molfile) != 0:
                    cpd = compound.Compound.create(molfile)
                    compound_dict[cpd.compound_name] = cpd
            # add the kat to kegg compound
            if database_name == "KEGG":
                kcf_files = glob.glob(working_directory + "/sources/KEGG/kcf/*")
                for kcf_file in kcf_files:
                    if os.path.getsize(kcf_file) != 0:
                        kcf_compound = parser.create_compound_kcf(kcf_file)
                        if "cpd:" + kcf_compound.cpd.compound_name in compound_dict:
                            parser.add_kat(compound_dict["cpd:" + kcf_compound.cpd.compound_name], kcf_compound)

            reaction_list = []
            # construct reactions
            if database_name == "KEGG":
                reaction_directory = working_directory + "/sources/KEGG/reaction/"
                rclass_directory = working_directory + "/sources/KEGG/rclass/"
                if not os.path.exists(reaction_directory):
                    raise OSError("The directory {0} does not exist.".format(reaction_directory))
                if not os.path.exists(rclass_directory):
                    raise OSError("The directory {0} does not exist.".format(rclass_directory))
                atom_mappings = parser.create_atom_mappings(rclass_directory, compound_dict)
                reaction_list = parser.create_reactions(reaction_directory, compound_dict, atom_mappings)

            elif database_name == "MetaCyc":
                reaction_file = working_directory + "/sources/MetaCyc/reaction.dat"
                atom_mapping_file = working_directory + "/sources/MetaCyc/atom-mapping.dat"
                if not os.path.exists(reaction_file):
                    raise OSError("The file {0} does not exist.".format(reaction_file))
                if not os.path.exists(atom_mapping_file):
                    raise OSError("The file {0} does not exist.".format(atom_mapping_file))
                reaction_list = parser.create_reactions(reaction_file, atom_mapping_file, compound_dict)

                # update compound
            for cpd_name in compound_dict:
                cpd = compound_dict[cpd_name]
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
            compound_dict.append(tools.open_jsonpickle(from_directory + "compounds.json"))
            reaction_list.append(tools.open_jsonpickle(from_directory + "reactions.json"))

        compound_harmonization_manager = harmonization.harmonize_compound_list(compound_dict)
        # while we do reaction harmonization, we need to pay attention to compound harmonization without same structural representations.
        # this includes: R group, linear-circular-transformation, resonance.
        reaction_harmonization_manager = harmonization.harmonize_reaction_list(reaction_list, compound_harmonization_manager)
        tools.save_to_jsonpickle(reaction_harmonization_manager, working_directory + "/harmonized/{0}.json".format("_".join(database_names)))

    elif args['test1']:

        aromatic_manager = aromatics.AromaticManager.decode(
            tools.open_jsonpickle("/mlab/data/hji236/projects/MDH_test/kcf_aromatic_manager.json"))
        kegg_cpd_path = "/mlab/data/hji236/projects/MDH_test/standardized/KEGG/molfile_test1"
        meta_cpd_path = "/mlab/data/hji236/projects/MDH_test/standardized/MetaCyc/molfile_test"
        to_file = "/mlab/data/hji236/projects/MDH_test/harmonized_edge_list_1.json"
        kegg_molfiles = glob.glob(kegg_cpd_path + "/*")
        #meta_molfiles = glob.glob(meta_cpd_path + "/*")
        kegg_dict = function_multiprocess(kegg_molfiles, construct_cpd)
        print("kegg parsed")
        #meta_dict = function_multiprocess(meta_molfiles, construct_cpd)
        #print("metacyc parsed")
        #
        aromatic_detection = aromatic_manager.detect_aromatic_substructures
        print(aromatic_detection)
        with multiprocessing.Pool() as pool:
            results = pool.map(aromatic_detection, list(kegg_dict.values()) )
        #print(len(aromatic_manager.aromatic_substructures))
        #for cpd_name in kegg_dict:
        #    cpd = kegg_dict[cpd_name]
        #    aromatic_manager.detect_aromatic_substructures(cpd)
            #cpd.define_bond_stereochemistry()
            #cpd.curate_invalid_n()
        #for cpd_name in meta_dict:
        #    cpd = meta_dict[cpd_name]
        #    aromatic_manager.detect_aromatic_substructures(cpd)
        #    cpd.define_bond_stereochemistry()
        #    cpd.curate_invalid_n()
        #compound_harmonization_manager = harmonization.harmonize_compound_list([kegg_dict, meta_dict])
        #edge_list = compound_harmonization_manager.get_edge_list()
        #tools.save_to_json(edge_list, to_file)


    elif args['test']:

        # aromatic_manager = aromatics.AromaticManager.decode(tools.open_jsonpickle("/mlab/data/hji236/projects/MDH_test/aromatic_manager.json"))
        cpds = tools.open_jsonpickle("/mlab/data/hji236/projects/MDH_test/KEGG_kcf_aromatics.json")
        update_cpds = []
        for cpd in cpds:
            update_cpds.append((cpd["compound_name"], cpd["atoms"], cpd["bonds"]))
        aromatic_manager = aromatics.AromaticManager.decode(update_cpds)
        save_file = "/mlab/data/hji236/projects/MDH_test/kcf_aromatics.json"
        tools.save_to_jsonpickle(aromatic_manager.encode(), save_file)


        cpd_path = "/mlab/data/hji236/projects/MDH_test/standardized/KEGG/molfile"
        kcf_path = "/mlab/data/hji236/projects/MDH_test/sources/KEGG/kcf"
        rclass_path = "/mlab/data/hji236/projects/MDH_test/sources/KEGG/rclass/"
        to_file = "/mlab/data/hji236/projects/MDH_test/kegg_atom_mappings.json"
        # parser = parser_dict['KEGG']
        #
        # molfiles = glob.glob(cpd_path + "/*")
        # compound_dict = function_multiprocess(molfiles, construct_cpd)
        #
        # print("mofile construct")
        
        # kcf_files = glob.glob(kcf_path + "/*")
        # compound_dict_kcf = function_multiprocess(kcf_files, construct_kcf)
        # print("kcf construct")
        #
        # for cpd in compound_dict:
        #     if cpd[4:] in compound_dict_kcf:
        #         try:
        #             parser.add_kat(compound_dict[cpd], compound_dict_kcf[cpd[4:]])
        #         except Exception as e:
        #             print(e)
        #             print("modified molfile, ", cpd)
        #             pass

                
        #       
        # for kcf_file in kcf_files:
        #     if os.path.getsize(kcf_file) != 0:
        #         kcf_compound = parser.create_compound_kcf(kcf_file)
        #         if "cpd:" + kcf_compound.cpd.compound_name in compound_dict:
        #             parser.add_kat(compound_dict["cpd:" + kcf_compound.cpd.compound_name], kcf_compound)

        # atom_mappings = parser.create_atom_mappings(rclass_path, compound_dict)
        # tools.save_to_json(atom_mappings, to_file)
        #
        # for cpd_name in compound_dict:
        #     cpd = compound_dict[cpd_name]
        #     #aromatic_manager.detect_aromatic_substructures(cpd)
        #     cpd.define_bond_stereochemistry()
        #     cpd.curate_invalid_n()



        # file1 = "/mlab/data/hji236/projects/MDH_test/KEGG_kcf_aromatics.json"
        # file2 = "/mlab/data/hji236/projects/MDH_test/indigo_aromatics.json"
        # combined = "/mlab/data/hji236/projects/MDH_test/combined_aromatics.json"
        # manager1 = aromatics.AromaticManager.decode(tools.open_jsonpickle(file1))
        # manager2 = aromatics.AromaticManager.decode(tools.open_jsonpickle(file2))
        # manager1.add_aromatic_substructures(manager2.aromatic_substructures)
        # print(len(manager1.aromatic_substructures))
        # tools.save_to_jsonpickle(manager1.encode(), combined)

        # file = "/mlab/data/hji236/projects/MDH_test/standardized/cpd:C00032.mol"
        # data = tools.open_text(file)
        # # print(data)
        # cpd = compound.Compound.create(file)
        # new_cpd = compound.Compound(cpd.compound_name, [atom.clone() for atom in cpd.atoms], [bond.clone() for bond in cpd.bonds])
        # atom_cnt = 0
        # for atom in cpd.atoms:
        #     if atom.in_cycle:
        #         atom_cnt += 1
        #         print(atom.default_symbol, atom.atom_number)
        # print(atom_cnt)
        # bond_cnt = 0
        # for bond in cpd.bonds:
        #     if cpd.atoms[bond.first_atom_number].in_cycle and cpd.atoms[bond.second_atom_number].in_cycle:
        #         bond_cnt += 1
        #         print(bond.first_atom_number, bond.second_atom_number, bond.bond_type)
        # print(bond_cnt)
        # file = "/mlab/data/hji236/projects/MDH_test/sources/KEGG/kcf/cpd:C00032"
        # parser = parser_dict['KEGG']
        # aromatic_manager = aromatics.AromaticManager()
        # kcf_cpd = parser.create_compound_kcf(file)
        # aromatic_manager.kegg_aromatize(kcf_cpd)
        # cpd = aromatic_manager.aromatic_substructures[0]
        # for atom in cpd.atoms:
        #     atom.in_cycle = False
        # new_cpd = compound.Compound(cpd.compound_name, cpd.atoms, cpd.bonds)
        # encoded = aromatic_manager.encode()
        # decoded = aromatics.AromaticManager.decode(encoded)




