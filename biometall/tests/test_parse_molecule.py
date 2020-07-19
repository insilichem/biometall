import pytest
from pytest import approx
from pathlib import Path
import numpy as np
import os
from biometall.modules import pdb

TEST_DATA_DIR = os.path.join(Path(__file__).resolve().parent, 'data')

testdata = [
            ('1oi0.pdb', '.pdb', [-1.1658,39.4787,-46.2727], 29.3750, 
                108, [4.701,58.174,-45.122],
                108, [0,0,0],
                108, [5.414,57.399,-44.034],
                108, [5.422,58.076,-46.424],
                108, [6.498,56.855,-44.259],
                108, 0, '1:A',
                'GLY', {'C', 'CA', 'N', 'O'}),
]

@pytest.mark.parametrize("file_name,extension,centroid,max_distance, \
                            len_alpha,first_alpha, \
                            len_beta,first_beta, \
                            len_carbon,first_carbon, \
                            len_nitrogen,first_nitrogen, \
                            len_oxygen,first_oxygen, \
                            len_res,column,res, \
                            name_for_res,atoms_in_res", testdata)

def test_parse_molecule(file_name, extension, centroid, max_distance, 
                        len_alpha, first_alpha,
                        len_beta, first_beta,
                        len_carbon, first_carbon,
                        len_nitrogen, first_nitrogen,
                        len_oxygen, first_oxygen,
                        len_res, column, res, 
                        name_for_res, atoms_in_res):
    file_path = os.path.join(TEST_DATA_DIR, file_name)

    with open(file_path, "r") as f:
        lines = f.read().splitlines()

    result = pdb._parse_molecule(lines, extension)
    assert centroid == approx(np.array(result[0]), abs=1e-3)
    assert max_distance == approx(result[1], abs=1e-3)
    assert len_alpha == len(result[2])
    assert first_alpha == approx(np.array(result[2][0]), abs=1e-3)
    assert len_beta == len(result[3])
    assert first_beta == approx(np.array(result[3][0]), abs=1e-3)
    assert len_carbon == len(result[4])
    assert first_carbon == approx(np.array(result[4][0]), abs=1e-3)
    assert len_nitrogen == len(result[5])
    assert first_nitrogen == approx(np.array(result[5][0]), abs=1e-3)
    assert len_oxygen == len(result[6])
    assert first_oxygen == approx(np.array(result[6][0]), abs=1e-3)
    assert len_res == len(result[7])
    assert len_res == len(result[8])
    assert column == result[7][res]
    assert res == result[8][column]
    assert name_for_res == result[9][res]
    assert atoms_in_res == result[10][res]

test_exception = [('1oi0_exception.pdb', '.pdb')]

@pytest.mark.parametrize("file_name,extension", test_exception)

def test_parse_molecule_exception(file_name, extension):
    file_path = os.path.join(TEST_DATA_DIR, file_name)

    with open(file_path, "r") as f:
        lines = f.read().splitlines()

    with pytest.raises(Exception):
        pdb._parse_molecule(lines, extension)