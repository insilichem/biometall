"""
pdb.py
Module with functions to manage the parsing and generation of .pdb files/structures
"""

import numpy as np
from .data import CONVERT_RES_NAMES, DIST_PROBE_ALPHA

def _parse_molecule(lines, file_extension):
    """
    Parses a protein and generates the necessary data for the calculation.

    All the atoms of the protein areparsed by searching those of type alpha
    carbon, beta carbon, backbone carbon, backbone nitrogen and backbone oxygen.
    Their coordinates are stored in separated numpy arrays for further use in
    the BioMetAll calculation. Also, it generates dictionaries to control the
    correspondence between residues and column numbers of the numpy arrays.
    Finally, the centroid and distance to the furthest atom of the protein are
    calculated to allow the generation of the grid of probes.

    Parameters
    ----------
    lines : array_like
        Array of str containing the structure of the protein
    file_extension : str
        Extension indicating the format of the `lines`

    Returns
    -------
    np.array
        3-float numpy array containing the centroid of the protein
    float
        Distance to the furthest atom, adding a security margin to account for
        a possible superficial coordination
    np.array
        Array of 3-D coordinates for all the alpha carbons of the protein
    np.array
        Array of 3-D coordinates for all the beta carbons of the protein
    np.array
        Array of 3-D coordinates for all the backbone carbons of the protein
    np.array
        Array of 3-D coordinates for all the backbone nitrogens of the protein
    np.array
        Array of 3-D coordinates for all the backbone oxygens of the protein
    np.array
        Array of 3-D coordinates for all the side-chain atoms of the protein
    dict
        Correspondence between number:chain of residue and column number
    dict
        Correspondence between column number and number of residue:chain
    dict
        Name of the residue given its number:chain
    dict
        Name of atoms contained in a given residue (indexed by number_res:chain)
    """
    if file_extension == '.pdb':
        #Extract residue information and assign column
        i = 0
        column_for_res = {}
        res_for_column = {}
        name_for_res = {}
        atoms_in_res = {}
        for line in lines:
            record_type = line[0:6]
            if record_type == "ATOM  ":
                atom_fullname = line[12:16]
                # get rid of whitespace in atom names
                split_list = atom_fullname.split()
                if len(split_list) != 1:
                    # atom name has internal spaces, e.g. " N B ", so
                    # we do not strip spaces
                    atom_name = atom_fullname
                else:
                    # atom name is like " CA ", so we can strip spaces
                    atom_name = split_list[0]

                if atom_name in ['CA', 'CB', 'C', 'N', 'O']:
                    altloc = line[16]
                    chainid = line[21]
                    resid = line[22:26].split()[0]
                    res = str(resid) + ":" + str(chainid)
                    resname = line[17:20]
                    if resname in list(CONVERT_RES_NAMES):
                        resname = CONVERT_RES_NAMES[resname]
                    if res not in list(column_for_res):
                        column_for_res[res] = i
                        res_for_column[i] = res
                        name_for_res[res] = resname
                        atoms_in_res[res] = set()
                        i += 1
                    atoms_in_res[res].add(atom_name)

        #Extract coordinates and atoms information
        alphas = [[0.0, 0.0, 0.0] for i in range(0, len(list(column_for_res)))]
        betas = [[0.0, 0.0, 0.0] for i in range(0, len(list(column_for_res)))]
        carbons = [[0.0, 0.0, 0.0] for i in range(0, len(list(column_for_res)))]
        nitrogens = [[0.0, 0.0, 0.0] for i in range(0, len(list(column_for_res)))]
        oxygens = [[0.0, 0.0, 0.0] for i in range(0, len(list(column_for_res)))]
        side_chains = []
        coords_array = [] #For calculate grid size

        for line in lines:
            record_type = line[0:6]
            if record_type == "ATOM  ":
                atom_fullname = line[12:16]
                # get rid of whitespace in atom names
                split_list = atom_fullname.split()
                if len(split_list) != 1:
                    # atom name has internal spaces, e.g. " N B ", so
                    # we do not strip spaces
                    atom_name = atom_fullname
                else:
                    # atom name is like " CA ", so we can strip spaces
                    atom_name = split_list[0]

                chainid = line[21]
                resid = line[22:26].split()[0]
                res = str(resid) + ":" + str(chainid)

                # atomic coordinates
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except Exception:
                    raise Exception("Invalid or missing coordinate(s) at \
                                    residue %s, atom %s" % (res, name))
                coord = [x, y, z]

                if atom_name == "CA":
                    # Coordinates for the grid
                    coords_array.append(coord)
                    # Coordinates for searching sites
                    alphas[column_for_res[res]] = coord
                elif atom_name == "CB":
                    # Coordinates for searching sites
                    betas[column_for_res[res]] = coord
                elif atom_name == "C":
                    # Coordinates for searching sites
                    carbons[column_for_res[res]] = coord
                elif atom_name == "N":
                    # Coordinates for searching sites
                    nitrogens[column_for_res[res]] = coord
                elif atom_name == "O":
                    # Coordinates for searching sites
                    oxygens[column_for_res[res]] = coord
                else: # Atom belongs to a side-chain
                    # Coordinates for discarding clashes
                    side_chains.append(coord)

        coords_array = np.array(coords_array)
        centroid =  np.mean(coords_array, axis=0)
        max_distance  = np.max(np.linalg.norm(coords_array - centroid, axis=1)) \
                        + DIST_PROBE_ALPHA['ALL'][1]

        alphas = np.array(alphas)
        betas = np.array(betas)
        carbons = np.array(carbons)
        nitrogens = np.array(nitrogens)
        oxygens = np.array(oxygens)
        side_chains = np.array(side_chains)
    return centroid, max_distance, alphas, betas, carbons, nitrogens, \
            oxygens, column_for_res, res_for_column, name_for_res, \
            atoms_in_res, side_chains

def _print_pdb(sorted_data, filename):
    """
    Generates a .pdb file containing the probes of the BioMetAll calculation.

    Each coordinating environment obtained from the BioMetAll calculation is
    stored as a different residue, being the centroid probe a Helium atom and
    all the probes Xeon atoms.

    Parameters
    ----------
    sorted_data : array_like
        Data obtained from the clustering of BioMetAll results
    filename : str
        Name of the output .pdb file. Usually with format `probes_xxxx.pdb` in
        the working directory.
    """
    file_pdb = open(filename,"w")
    num_at = 0
    num_res = 0
    for one_result in sorted_data:
        chains = set()
        for r in one_result[0]:
            r = r.strip("_BCK")
            chains.add(r.split(":")[1])
        cen_str = ""
        for r in one_result[1]:
            crd_center = "{:.8s}".format(str(round(float(r),3)))
            if len(crd_center)<8:
                crd_center = " "*(8-len(crd_center)) + crd_center
                cen_str += crd_center
            else:
                cen_str += crd_center
        num_at += 1
        num_res += 1
        for ch in chains:
            file_pdb.write("ATOM" +" "*(7-len(str(num_at))) + "%s  HE  SLN %s" %(num_at, ch))
            file_pdb.write(" "*(3-len(str(num_res))) + "%s     %s  1.00  0.00          HE\n" %(num_res, cen_str))
        for prob in one_result[4]:
            num_at += 1
            prb_str = ""
            for p in prob:
                prb_center = "{:.8s}".format(str(round(float(p),3)))
                if len(prb_center)<8:
                    prb_center = " "*(8-len(prb_center)) + prb_center
                    prb_str += prb_center
                else:
                    prb_str += prb_center
            for ch in chains:
                file_pdb.write("ATOM" +" "*(7-len(str(num_at))) + "%s  XE  SLN %s" %(num_at, ch))
                file_pdb.write(" "*(3-len(str(num_res))) + "%s     %s  1.00  0.00          XE\n" %(num_res, prb_str))
    file_pdb.close()
