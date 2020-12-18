"""
biometall.py
Software to predict feasible metal binding areas into proteins

Handles the primary functions
"""

import numpy as np
import sys
import time
import math
import os
import re
import multiprocessing
import psutil
import itertools
import copy
import urllib.request
from functools import partial
from ._version import get_versions

# BioMetAll modules
from .modules.data import DIST_PROBE_ALPHA, DIST_PROBE_BETA, ANGLE_PAB
from .modules.data import DIST_PROBE_OXYGEN, ANGLE_POC
from .modules.data import ALLOWED_FILE_TYPES
from .modules.grid import _chunk_size, _calculate_center_and_radius, _grid
from .modules.pdb import _parse_molecule, _print_pdb
from .modules.motif import _check_actual_motif, _check_possible_mutations

def run(inputfile, min_coordinators=3, min_sidechain=2,
        residues='[ASP,HIS,GLU,CYS]', motif='', grid_step=1.0,
        consider_backbone_residues='[]', cluster_cutoff=0.0, pdb=False,
        propose_mutations_to='', custom_radius=None, custom_center=None,
        cores_number=None, backbone_clashes_threshold=1.0,
        sidechain_clashes_threshold=0.0, cmd_str="", **kwargs):
    # Print header
    versions = get_versions()
    __version__ = versions['version']
    str_header = "You are running BioMetAll {}".format(__version__)
    print(str_header)
    print("*" * len(str_header))
    print("")

    filename, file_extension = os.path.splitext(inputfile)

    if motif:
        motifs = list(map(str, motif.strip('[]').split(',')))
        motifs_list = []
        for mot in motifs:
            motifs_list.append(list(map(str, mot.split('/'))))
        motifs = motifs_list
        motifs.sort(key=len)
    else:
        motifs = None
    if propose_mutations_to:
        mutated_motif = list(map(str,
                                propose_mutations_to.strip('[]').split(',')))
        mutated_motif_list = []
        for mot in mutated_motif:
            mutated_motif_list.append(list(map(str, mot.split('/'))))
        mutated_motif = mutated_motif_list
        mutated_motif.sort(key=len)
    else:
        mutated_motif = None

    #Set residues to consider as coordinating
    if propose_mutations_to:
        tot = mutated_motif_list + motifs_list
        residues = list(set(x for l in tot for x in l))
    elif motif:
        tot = motifs_list
        residues = list(set(x for l in tot for x in l))
    else:
        residues = list(map(str, residues.strip('[]').split(',')))

    # Control of a correct usage of parameters
    if propose_mutations_to and motif:
        #min_coordinators should be at least len(motif) + len(mutated_motif)
        if min_coordinators < (len(motifs) + len(mutated_motif)):
            min_coordinators = len(motifs) + len(mutated_motif)
    elif propose_mutations_to and not motif:
        print("To propose mutations is necessary to set a base motif with --motif parameter.")
        sys.exit()
    elif motif:
        #min_coordinators should be at least len(motif)
        if min_coordinators < len(motifs):
            min_coordinators = len(motifs)
            print("min_coordinators has been set to {} due to the motif length".format(len(motifs)))

    #Set residues to consider as coordinating in backbone oxygens
    if consider_backbone_residues == '[]':
        consider_backbone_residues = None
        min_sidechain = min_coordinators
    elif consider_backbone_residues == 'ALL':
        consider_backbone_residues = 'ALL'
    else:
        consider_backbone_residues = list(map(str, consider_backbone_residues.strip('[]').split(',')))


    t0 = time.time()

    if not cores_number:
        cores_number = psutil.cpu_count(logical=False)
    pool = multiprocessing.Pool(cores_number)

    #Parse inputfile to obtain atom coords
    if file_extension in ALLOWED_FILE_TYPES:
        with open(inputfile, "r") as f:
            lines = f.read().splitlines()
    elif file_extension == "":
        try:
            response = urllib.request.urlretrieve('http://files.rcsb.org/download/%s.pdb' % (filename), '%s.pdb' % (filename))
            file_extension = ".pdb"
            with open(filename + file_extension, "r") as f:
                lines = f.read().splitlines()
        except:
            raise Exception("The input should be a valid PDB code or a file of a valid type: " + str(ALLOWED_FILE_TYPES))
    else:
        raise Exception("The input should be a valid PDB code or a file of a valid type: " + str(ALLOWED_FILE_TYPES))

    centroid, radius, alphas, betas, carbons, nitrogens, oxygens, column_for_res, res_for_column, name_for_res, atoms_in_res, side_chains = _parse_molecule(lines, file_extension)

    #Alpha-Beta and Oxygen-Carbon distances for further use
    alpha_beta_distances = np.sqrt((np.square(betas-alphas).sum(axis=1)))
    oxygen_carbon_distances = np.sqrt((np.square(carbons-oxygens).sum(axis=1)))

    #Construct grid
    if custom_center and custom_radius:
        centroid = np.array(list(map(float, custom_center.strip('[]').split(','))))
        radius = custom_radius
    grid = _grid(centroid, radius, grid_step)

    #Split grid in chunks (to ensure a good use of the memory/processors)
    n = _chunk_size(len(grid), len(alphas), cores_number)

    grids = [grid[i * n:(i + 1) * n] for i in range((len(grid) + n - 1) // n )]

    coordination_chunks = pool.map(partial(_test_chunk, alphas=alphas,
                        betas=betas, carbons=carbons, nitrogens=nitrogens,
                        oxygens=oxygens, side_chains=side_chains,
                        alpha_beta_distances=alpha_beta_distances,
                        oxygen_carbon_distances=oxygen_carbon_distances,
                        name_for_res=name_for_res,
                        column_for_res=column_for_res,
                        res_for_column=res_for_column,
                        atoms_in_res=atoms_in_res,
                        consider_backbone_residues=consider_backbone_residues,
                        DIST_PROBE_ALPHA=DIST_PROBE_ALPHA,
                        DIST_PROBE_BETA=DIST_PROBE_BETA,
                        ANGLE_PAB=ANGLE_PAB,
                        bck_clashes=backbone_clashes_threshold,
                        sc_clashes=sidechain_clashes_threshold), grids)

    centers, mutations = clustering(coordination_chunks, residues, motifs, min_coordinators, min_sidechain,
                    consider_backbone_residues, mutated_motif, cluster_cutoff, filename, name_for_res, res_for_column, atoms_in_res)
    #Print results
    sorted_data = sorted(centers, key=lambda x: x[2], reverse=True)
    if motif:
        file_name_addendum = "_" + motif.replace("[", "").replace("]", "").replace(",", "_").replace("/", "-")
    else:
        file_name_addendum = ""
    if pdb:
        pdb_filename = "probes_%s%s.pdb" %(os.path.basename(filename), file_name_addendum)
        pdb_filename = os.path.join(os.path.dirname(inputfile), pdb_filename)
        _print_pdb(sorted_data, pdb_filename)

    mutations_width, radius_width, coord_width, probes_width, residues_width, pos_width = len('Proposed mutations'), len('Radius search'), len('Coordinates of center'), len('Num. probes'), len('Coordinating residues'), 1
    lines=[]
    for pos, one_center in enumerate(sorted_data,1):
        residues = []
        for res in one_center[0]:
            if not "_BCK" in res:
                residues.append(str(name_for_res[res]) + ":" + res)
            else:
                res_bck = res.split("_")[0]
                residues.append(str(name_for_res[res_bck]) + ":" + res_bck + "_BCK")
                continue

        residues_str = ' '.join(res for res in residues)
        coord_str = ' '.join([str(format(r, '.3f')) for r in one_center[1]])
        pos_str, radius_str, probes_str = str(pos), str(format(one_center[3], '.3f')), str(one_center[2])

        if propose_mutations_to:
            mutations_str = ' '.join([str(res) + ":" + str(mut) for (res,mut) in mutations[one_center[0]]])
            if len(mutations_str) > mutations_width:
                mutations_width = len(mutations_str)
        else:
            mutations_str = ''
        if len(pos_str) > pos_width:
            pos_width = len(pos_str)
        if len(coord_str) > coord_width:
            coord_width = len(coord_str)
        if len(probes_str) > probes_width:
            probes_width = len(probes_str)
        if len(residues_str) > residues_width:
            residues_width = len(residues_str)
        if len(radius_str) > radius_width:
            radius_width = len(radius_str)
        lines.append((pos_str, residues_str, coord_str, probes_str, radius_str, mutations_str))
    print(' {:>{pos_width}} | {:^{residues_width}} | {:{coord_width}} | {:{probes_width}} | {:{radius_width}} | {:{mutations_width}}'.format(
            '#', 'Coordinating residues', 'Coordinates of center', 'Num. probes', 'Radius search' , 'Proposed mutations', pos_width=pos_width, residues_width=residues_width, coord_width=coord_width,
            probes_width=probes_width, radius_width=radius_width, mutations_width=mutations_width))
    print('-{}-+-{}-+-{}-+-{}-+-{}-+-{}-'.format('-'*pos_width, '-'*residues_width, '-'*coord_width, '-'*probes_width, '-'*radius_width, '-'*mutations_width))
    for line in lines:
        print(' {:>{pos_width}} | {:<{residues_width}} | {:^{coord_width}} | {:>{probes_width}} | {:<{radius_width}} | {:<{mutations_width}} '.format(
                line[0], line[1], line[2], line[3], line[4], line[5], pos_width=pos_width, residues_width=residues_width,
                coord_width=coord_width, probes_width=probes_width, radius_width=radius_width, mutations_width=mutations_width))

    text_filename = "results_biometall_%s%s.txt" %(os.path.basename(filename), file_name_addendum)
    text_filename = os.path.join(os.path.dirname(inputfile), text_filename)
    f = open(text_filename, "w")
    str_header = "*****BioMetAll {}".format(__version__)
    f.write(str_header)
    f.write('\n')
    f.write(cmd_str)
    f.write('\n')
    f.write(' {:>{pos_width}} | {:^{residues_width}} | {:{coord_width}} | {:{probes_width}} | {:{radius_width}} | {:{mutations_width}} \n'.format(
            '#', 'Coordinating residues', 'Coordinates of center', 'Num. probes', 'Radius search' , 'Proposed mutations', pos_width=pos_width, residues_width=residues_width, coord_width=coord_width,
            probes_width=probes_width, radius_width=radius_width, mutations_width=mutations_width))
    f.write('-{}-+-{}-+-{}-+-{}-+-{}-+-{}-\n'.format('-'*pos_width, '-'*residues_width, '-'*coord_width, '-'*probes_width, '-'*radius_width, '-'*mutations_width ))
    for line in lines:
        f.write(' {:>{pos_width}} | {:<{residues_width}} | {:^{coord_width}} | {:>{probes_width}} | {:<{radius_width}} | {:<{mutations_width}} \n'.format(
                line[0], line[1], line[2], line[3], line[4], line[5], pos_width=pos_width, residues_width=residues_width,
                coord_width=coord_width, probes_width=probes_width, radius_width=radius_width, mutations_width=mutations_width))
    f.write("*****Calculation took {0:.2f} seconds".format(time.time() - t0))
    f.close()
    print("{0:.2f} seconds".format(time.time() - t0))

def clustering(coordination_chunks, residues, motifs, min_coordinators, min_sidechain,
            consider_backbone_residues, mutated_motif, cluster_cutoff, filename,
            name_for_res, res_for_column, atoms_in_res):
    dict_cluster = {}
    dict_mutations = {}
    centers = []
    if motifs:
        motif_possibilities = list(itertools.product(*motifs))
    if mutated_motif:
        mutated_motif_possibilities = list(itertools.product(*mutated_motif))
        residues_of_mutated_motif = set(x for l in mutated_motif for x in l)

    for probes, coordinations, discarded in coordination_chunks:
        coordinators = {}
        for possible_coord_name in residues: #Sidechain coordinations
            for probe_idx, res_idx in coordinations[possible_coord_name][0,:,:]:
                if not mutated_motif and (possible_coord_name != name_for_res[res_for_column[res_idx]]): #discarded due to different residue name
                    continue
                if probe_idx in discarded: #residue is discarded for that probe
                    continue
                elif ('CA' not in atoms_in_res[res_for_column[res_idx]]) or ('CB' not in atoms_in_res[res_for_column[res_idx]]):
                    continue
                else: #coordination is valid for that probe and residue
                    if probe_idx not in list(coordinators):
                        coordinators[probe_idx] = {res_idx: [possible_coord_name]}
                    else:
                        if res_idx not in list(coordinators[probe_idx]):
                            coordinators[probe_idx][res_idx] = [possible_coord_name]
                        else:
                            coordinators[probe_idx][res_idx].append(possible_coord_name)
        if consider_backbone_residues:
            for probe_idx, res_idx in coordinations["BCK"][0,:,:]:
                if probe_idx in discarded:
                    continue
                if (name_for_res[res_for_column[res_idx]] not in consider_backbone_residues) and consider_backbone_residues != "ALL": #discarded due to not searched residue
                    continue
                else: #coordination is valid for that probe and residue
                    if probe_idx not in list(coordinators):
                        coordinators[probe_idx] = {res_idx: [name_for_res[res_for_column[res_idx]] + "_BCK"]}
                    else:
                        if res_idx not in list(coordinators[probe_idx]):
                            coordinators[probe_idx][res_idx] = [name_for_res[res_for_column[res_idx]] + "_BCK"]
                        else:
                            coordinators[probe_idx][res_idx].append(name_for_res[res_for_column[res_idx]] + "_BCK")
        for probe_idx in list(coordinators):
            if motifs:
                #A minimum motif should be accomplished without mutations (contained in motifs/motif_possibilities)
                actual_motif_solutions = _check_actual_motif(motif_possibilities, coordinators[probe_idx], name_for_res, res_for_column)
            if mutated_motif:
                #The resting coordinators of each solution should be mutable to complete the propose_mutations_to motif
                mutation_motif_solutions = {}
                for actual_motif_solution in actual_motif_solutions:
                    resting_coordinators = copy.deepcopy(coordinators[probe_idx])
                    for r in actual_motif_solution:
                        del resting_coordinators[r]
                    mutation_motif_solutions = _check_possible_mutations(mutated_motif_possibilities, residues_of_mutated_motif, resting_coordinators)

            if motifs and mutated_motif and actual_motif_solutions and mutation_motif_solutions:
                #Searching for already present motifs plus proposing mutations
                #already present motifs
                for actual_motif_solution in actual_motif_solutions:
                    item = [probes[probe_idx], tuple(sorted([res_for_column[c] for c in actual_motif_solution]))]
                    if item[1] not in list(dict_cluster):
                        dict_cluster[item[1]] = [item[0]]
                    else:
                        dict_cluster[item[1]] += [item[0]]
                    #Mutation information
                    for m in list(mutation_motif_solutions):
                        if item[1] not in dict_mutations:
                            dict_mutations[item[1]] = {res_for_column[m]: [[el, 1] for el in mutation_motif_solutions[m]]}
                        else:
                            if res_for_column[m] not in dict_mutations[item[1]]:
                                dict_mutations[item[1]][res_for_column[m]] = [[el, 1] for el in mutation_motif_solutions[m]]
                            else:
                                for el1 in mutation_motif_solutions[m]:
                                    if el1 not in [el2[0] for el2 in dict_mutations[item[1]][res_for_column[m]]]:
                                        dict_mutations[item[1]][res_for_column[m]].append([el1, 1])
                                    else:
                                        for el2 in dict_mutations[item[1]][res_for_column[m]]:
                                            if el1 == el2[0]:
                                                el2[1] += 1
                                                break
            elif motifs and not mutated_motif:
                #Searching for already present motifs
                for actual_motif_solution in actual_motif_solutions:
                    item = [probes[probe_idx], tuple(sorted([res_for_column[c] for c in actual_motif_solution]))]
                    if item[1] not in list(dict_cluster):
                        dict_cluster[item[1]] = [item[0]]
                    else:
                        dict_cluster[item[1]] += [item[0]]
            elif not motifs and not mutated_motif:
                #Searching number of coordinators
                sidechain_coord = []
                bck_coord = []
                for item in list(coordinators[probe_idx]):
                    actual_residue_name = name_for_res[res_for_column[item]]
                    if actual_residue_name in coordinators[probe_idx][item]:
                        sidechain_coord.append(res_for_column[item])
                    if actual_residue_name + "_BCK" in coordinators[probe_idx][item]:
                        bck_coord.append(str(res_for_column[item]) + "_BCK")
                if ((len(sidechain_coord) + len(bck_coord)) >= min_coordinators) and (len(sidechain_coord) >= min_sidechain):
                    c = tuple(sorted(sidechain_coord + bck_coord))
                    for coord_environment in list(dict_cluster):
                        if set(c).issuperset(coord_environment):
                            dict_cluster[coord_environment] += [probes[probe_idx]]
                    coordination_possibilities = []
                    for i in range(min_coordinators, len(c)+1):
                        coordination_possibilities += list(itertools.combinations(c,i))
                    coordination_possibilities = [tuple(sorted(el)) for el in coordination_possibilities]
                    for el in coordination_possibilities:
                        sc_num = sum("BCK" not in L for L in el)
                        if (el not in dict_cluster) and sc_num >= min_sidechain:
                            dict_cluster[el] = [probes[probe_idx]]

    #Order Mutation information
    for coord_environment in list(dict_mutations):
        for residue in list(dict_mutations[coord_environment]):
            dict_mutations[coord_environment][residue].sort(key = lambda x: x[1], reverse=True)   #order every residue by number of probes
    for coord_environment in list(dict_mutations):
        dict_mutations[coord_environment] = sorted(dict_mutations[coord_environment].items(), key=lambda item: item[1][0][1], reverse=True)

    try:
        max_probes = max([len(v) for v in dict_cluster.values()])
    except:
        print("None possible coordinating sites have been found. Try again with other parameters or check/change the input file.")
        sys.exit()

    for coord_residues,probes in dict_cluster.items():
        if len(probes) >= max_probes*cluster_cutoff:
            center, radius_search = _calculate_center_and_radius(probes)
            centers.append((coord_residues, center, len(probes), radius_search, probes))

    return centers, dict_mutations

def _test_chunk(grid, alphas, betas, carbons, nitrogens, oxygens, side_chains,
                        alpha_beta_distances, oxygen_carbon_distances,
                        name_for_res, column_for_res, res_for_column,
                        atoms_in_res, consider_backbone_residues,
                        DIST_PROBE_ALPHA, DIST_PROBE_BETA, ANGLE_PAB,
                        bck_clashes, sc_clashes):
    alpha_distances = np.sqrt((np.square(grid[:,np.newaxis]-alphas).sum(axis=2)))
    beta_distances = np.sqrt((np.square(grid[:,np.newaxis]-betas).sum(axis=2)))
    carbon_distances = np.sqrt((np.square(grid[:,np.newaxis]-carbons).sum(axis=2)))
    nitrogen_distances = np.sqrt((np.square(grid[:,np.newaxis]-nitrogens).sum(axis=2)))
    oxygen_distances = np.sqrt((np.square(grid[:,np.newaxis]-oxygens).sum(axis=2)))

    PAB_angles = np.arccos((np.square(alpha_distances) + np.square(alpha_beta_distances) - np.square(beta_distances)) / (2*alpha_distances*alpha_beta_distances))
    POC_angles = np.arccos((np.square(oxygen_distances) + np.square(oxygen_carbon_distances) - np.square(carbon_distances)) / (2*oxygen_distances*oxygen_carbon_distances))

    coords = {}
    # include backbone coordinations
    if consider_backbone_residues:
        coords["BCK"] = np.dstack(np.where((DIST_PROBE_OXYGEN[0]<=oxygen_distances) & (oxygen_distances<=DIST_PROBE_OXYGEN[1]) &
                                (ANGLE_POC[0]<=POC_angles) & (POC_angles<=ANGLE_POC[1])))
    # include sidechain coordinations (Alpha+Beta)
    for res_name in list(DIST_PROBE_ALPHA):
        coords[res_name] = np.dstack(np.where((DIST_PROBE_ALPHA[res_name][0]<=alpha_distances) & (alpha_distances<=DIST_PROBE_ALPHA[res_name][1]) &
                                                (DIST_PROBE_BETA[res_name][0]<=beta_distances) & (beta_distances<=DIST_PROBE_BETA[res_name][1]) &
                                                (ANGLE_PAB[res_name][0]<=PAB_angles) & (PAB_angles<=ANGLE_PAB[res_name][1])))
    # If there is a clash (distance < bck_clashes) with a backbone atom,
    # no coordination is possible for that probe
    discarded = set(np.where((oxygen_distances<bck_clashes) |
                            (carbon_distances<bck_clashes) |
                            (nitrogen_distances<bck_clashes) |
                            (alpha_distances<bck_clashes))[0])

    if sc_clashes > 0:
        for sc_atom in side_chains:
            sidechain_distances = np.sqrt((np.square(grid[:,np.newaxis]-sc_atom).sum(axis=2)))
            discarded = discarded.union(set(np.where(sidechain_distances<sc_clashes)[0]))

    return [grid, coords, discarded]
