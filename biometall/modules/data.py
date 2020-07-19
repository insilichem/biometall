"""
data.py
Module to store the different values used in BioMetAll calculations
"""

#Dictionary with ranges for Metal-Alpha carbon distances (in angstroms)
DIST_PROBE_ALPHA = {
    'ASP': (3.907, 6.192),
    'HIS': (3.076, 7.098),
    'GLU': (3.591, 8.303),
    'CYS': (3.248, 5.451),
    'ASN': (4.048, 5.701),
    'THR': (3.313, 5.096),
    'SER': (3.384, 5.171),
    'GLN': (3.217, 7.835),
    'MET': (4.426, 6.706),
    'TYR': (6.896, 9.016),
    'LYS': (5.134, 9.280),
    'ARG': (3.409, 10.162),
    'LEU': (3.784, 4.678),
    'TRP': (3.912, 7.478),
    'ILE': (4.236, 4.817),
    'VAL': (3.517, 5.003),
    'ALA': (3.402, 5.126),
    'PHE': (3.499, 4.922),
    'PRO': (3.646, 4.529),
    'ALL': (3.390, 7.218)
}

#Dictionary with ranges for Metal-Beta carbon distances (in angstroms)
DIST_PROBE_BETA = {
    'ASP': (3.658, 5.052),
    'HIS': (3.047, 6.073),
    'GLU': (4.123, 6.108),
    'CYS': (2.794, 3.829),
    'ASN': (3.866, 5.250),
    'THR': (2.530, 4.371),
    'SER': (2.587, 4.353),
    'GLN': (3.512, 6.565),
    'MET': (3.654, 5.094),
    'TYR': (6.367, 8.087),
    'LYS': (4.358, 7.997),
    'ARG': (2.913, 9.042),
    'LEU': (4.353, 5.622),
    'TRP': (4.244, 6.573),
    'ILE': (5.039, 5.702),
    'VAL': (4.358, 5.722),
    'ALA': (4.254, 5.991),
    'PHE': (4.289, 5.590),
    'PRO': (4.860, 5.246),
    'ALL': (2.666, 6.351)
}

#Dictionary with ranges for Metal-Alpha carbon-Beta carbon angles (in Rad)
ANGLE_PAB = {
    'ASP': (0.003, 1.871),
    'HIS': (0.001, 2.018),
    'GLU': (0.000, 2.305),
    'CYS': (0.006, 1.734),
    'ASN': (0.682, 1.697),
    'THR': (0.184, 1.647),
    'SER': (0.163, 1.521),
    'GLN': (0.012, 2.268),
    'MET': (0.016, 1.228),
    'TYR': (0.501, 1.433),
    'LYS': (0.029, 1.918),
    'ARG': (0.024, 1.991),
    'LEU': (1.464, 2.470),
    'TRP': (0.471, 2.703),
    'ILE': (1.826, 2.169),
    'VAL': (1.639, 2.166),
    'ALA': (1.584, 2.368),
    'PHE': (1.666, 2.300),
    'PRO': (1.746, 2.678),
    'ALL': (0.092, 1.601)
}

# Range for distance Metal-Backbone Oxygen (in angstroms)
DIST_PROBE_OXYGEN = (1.809, 3.000)

# Range for angle Metal-Backbone oxygen-Backbone carbon (in Rad)
ANGLE_POC = (1.774, 3.139)

# Type of files allowed in BioMetAll calculations
ALLOWED_FILE_TYPES = ('.pdb'),

# Conversion of some particular amino acid names to standard ones
CONVERT_RES_NAMES = {
    'CYX': 'CYS',
    'CYM': 'CYS',
    'ASH': 'ASN',
    'GLH': 'GLN',
    'HIE': 'HIS',
    'HID': 'HIS',
    'HIP': 'HIS',
    'HYP': 'PRO',
    'LYN': 'LYS'
}