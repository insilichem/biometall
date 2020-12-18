import argparse
import multiprocessing
from biometall import run
import sys

def parse_cli():
    p = argparse.ArgumentParser()
    p.add_argument('inputfile', type=str,
        help='Molecule pdb file to be analysed')
    p.add_argument('--min_coordinators', type=int, default=3,
        help='Minimum number of coordinating residues of a given grid probe for that probe to be considered as potentially coordinating. Default: 3')
    p.add_argument('--residues', type=str, default='[ASP,HIS,GLU,CYS]',
        help='Residues that are potentially coordinating. Requires a list of 3-letter aa code separated by commas. Default: [ASP,CYS,GLU,HIS]')
    p.add_argument('--motif', type=str, default='',
        help='Search for concrete potentially coordinating motifs. Check the manual for syntax (e.g. [HIS,HIS,GLU/ASP]), Default: None')
    p.add_argument('--mutations', type=str, dest='propose_mutations_to', default='',
        help='In case of searching a motif, propose possible mutations to complete the motif. Check the manual for syntax (e.g. [ASP/GLU]). Default: None')
    p.add_argument('--grid', type=float, dest='grid_step', default=1.0,
        help='Distance between two probes of the grid defining the search zone. Default: 1.0')
    p.add_argument('--center', type=str, dest='custom_center', default='',
        help='Center of the search sphere. Only if you want to scan a restricted zone of the molecule. Default: None')
    p.add_argument('--radius', type=float, dest='custom_radius', default=None,
        help='Radius of the search sphere. Only if you want to scan a restricted zone of the molecule. Default: None')
    p.add_argument('--cutoff', type=float, dest='cluster_cutoff', default=0,
        help='Cluster density cutoff. The possible binding site will only be printed if the cluster has more probes than the cutoff*most populated cluster. It should be between 0 and 1.0. Default: 0')
    p.add_argument('--pdb', action='store_true', default=False,
        help='Output a probes.pdb file with all possible binding zones. Default: False')
    p.add_argument('--cores', type=int, dest='cores_number', default=0,
        help='Set how many cores to use (if not set, all physical cores of the CPU will be used). Default: all physical cores')
    p.add_argument('--backbone', type=str, dest='consider_backbone_residues', default='[]',
        help='Consider possible coordinations with backbone oxygens in the requested residues list. Requires a list of 3-letter aa codes separated by commas, or "ALL" to consider all residues as possible coordinators. Default: None')
    p.add_argument('--min_sidechain', type=int, default=2,
        help='Minimum number of sidechain coordinators that a solution could have. It should be used together with --backbone parameter to control the type of coordinations that appear in the solutions. Default: 2')
    p.add_argument('--backbone_clashes', type=float,
        dest='backbone_clashes_threshold', default=1.0,
        help='Distance from a grid probe to a backbone atom that defines a clash. The probes at less distance from any backbone atom will be discarded. For example, if set to 1.0, all probes nearer than 1.0 Angstroms from any backbone atom of the protein will be discarded. Default: 1.0')
    p.add_argument('--sidechain_clashes', type=float,
        dest='sidechain_clashes_threshold', default=0.0,
        help='Distance from a grid probe to a sidechain atom that defines a clash. The probes at less distance from any sidechain atom will be discarded. For example, if set to 1.0, all probes nearer than 1.0 Angstroms from any sidechain atom of the protein will be discarded. Default: 0.0')
    cmd_str = "*****biometall " + " ".join(sys.argv[1:])

    return p.parse_args(), cmd_str

def main():
    args, cmd_str = parse_cli()
    run(**vars(args), cmd_str=cmd_str)

if __name__ == '__main__':
    multiprocessing.freeze_support()
    main()
