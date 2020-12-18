import pytest
import tempfile
from shutil import copyfile
from pathlib import Path
import os

import biometall

TEST_DATA_DIR = os.path.join(Path(__file__).resolve().parent, 'data')

testdata = [
            ('1oi0.pdb', 3, 2, '[ASP,HIS,GLU,CYS]', '', 1.0, '[]', 0.0,
                True, '', None, None, None, 1.0, 0.0,
                '1oi0_output0.txt', 'results_biometall_1oi0.txt',
                '1oi0_output0.pdb', 'probes_1oi0.pdb'),
            ('1oi0.pdb', 3, 2, '[ASP,HIS,GLU,CYS]', '[GLU,GLU,ASP/TYR]', 1.0,
                '[]', 0.2, True, '', None, None, None, 1.0, 0.0,
                '1oi0_output1.txt', 'results_biometall_1oi0_GLU_GLU_ASP-TYR.txt',
                '1oi0_output1.pdb', 'probes_1oi0_GLU_GLU_ASP-TYR.pdb'),
            ('1oi0.pdb', 3, 2, '[ASP,HIS,GLU,CYS]', '', 1.0,
                '[HIS]', 0.0, True, '', None, None, None, 1.0, 0.0,
                '1oi0_output2.txt', 'results_biometall_1oi0.txt',
                '1oi0_output2.pdb', 'probes_1oi0.pdb'),
            ('1oi0.pdb', 3, 2, '[ASP,HIS,GLU,CYS]', '[GLU,GLU,ASP]', 1.0,
                '[]', 0.0, True, '[TYR,TYR]', None, None, None, 1.0, 0.0,
                '1oi0_output3.txt', 'results_biometall_1oi0_GLU_GLU_ASP.txt',
                '1oi0_output3.pdb', 'probes_1oi0_GLU_GLU_ASP.pdb'),
            ('1oi0.pdb', 3, 2, '[ASP,HIS,GLU,CYS]', '', 1.0, '[]', 0.0,
                True, '', None, None, None, 1.5, 1.5,
                '1oi0_output4.txt', 'results_biometall_1oi0.txt',
                '1oi0_output4.pdb', 'probes_1oi0.pdb'),
]

CONTENT = "content"

@pytest.mark.parametrize("inputfile,min_coordinators,min_sidechain,residues, \
                            motif,grid_step,consider_backbone_residues, \
                            cluster_cutoff,pdb,propose_mutations_to, \
                            custom_radius,custom_center,cores_number, \
                            backbone_clashes_threshold, \
                            sidechain_clashes_threshold, \
                            output_text,output_text_name,output_pdb, \
                            output_pdb_name", testdata)

def test_full_calculation(inputfile, min_coordinators, min_sidechain, residues,
                            motif, grid_step, consider_backbone_residues,
                            cluster_cutoff, pdb, propose_mutations_to,
                            custom_radius, custom_center, cores_number,
                            backbone_clashes_threshold,
                            sidechain_clashes_threshold,
                            output_text, output_text_name, output_pdb,
                            output_pdb_name,tmp_path):

    d = Path(os.path.join(tmp_path, "full"))
    d.mkdir()

    inputfile_path = os.path.join(TEST_DATA_DIR, inputfile)
    tmp_inputfile = os.path.join(d, inputfile)
    copyfile(inputfile_path, tmp_inputfile)

    biometall.run(tmp_inputfile, min_coordinators, min_sidechain,
                    residues, motif, grid_step, consider_backbone_residues,
                    cluster_cutoff, pdb, propose_mutations_to, custom_radius,
                    custom_center, cores_number, backbone_clashes_threshold,
                    sidechain_clashes_threshold)

    output_text_path = os.path.join(TEST_DATA_DIR, output_text)
    result_text_path = os.path.join(d, output_text_name)
    with open(output_text_path, "r") as output_file, \
         open(result_text_path, "r") as result_file:
        assert output_file.readlines()[0:-1] == result_file.readlines()[2:-1]

    if output_pdb:
        output_pdb_path = os.path.join(TEST_DATA_DIR, output_pdb)
        result_pdb_path = os.path.join(d, output_pdb_name)
        with open(output_pdb_path, "r") as output_pdb, \
             open(result_pdb_path, "r") as result_pdb:
            assert output_pdb.readlines() == result_pdb.readlines()
