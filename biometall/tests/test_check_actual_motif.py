import pytest

from biometall.modules import motif


testdata = [
            ([('HIS', 'ASP'), ('HIS', 'GLU')],
                {69: ['ASP'], 80: ['HIS'], 95: ['TYR']}, 
                {'83:A': 'ASP', '94:A': 'HIS', '109:A': 'TYR'}, 
                {69: '83:A', 80: '94:A', 95: '109:A'},
                {(69, 80)}),
            ([('HIS', 'ASP'), ('HIS', 'GLU')],
                {69: ['SER'], 80: ['HIS'], 95: ['TYR']}, 
                {'83:A': 'SER', '94:A': 'HIS', '109:A': 'TYR'}, 
                {69: '83:A', 80: '94:A', 95: '109:A'},
                set()),
            ([('SER', 'TYR'), ('HIS', 'GLU')],
                {69: ['SER'], 80: ['HIS'], 95: ['TYR']}, 
                {'83:A': 'SER', '94:A': 'HIS', '109:A': 'TYR'}, 
                {69: '83:A', 80: '94:A', 95: '109:A'},
                {(69, 95)}),
            ([('SER', 'TYR'), ('HIS', 'HIS')],
                {69: ['HIS'], 80: ['HIS'], 95: ['HIS']}, 
                {'83:A': 'HIS', '94:A': 'HIS', '109:A': 'HIS'}, 
                {69: '83:A', 80: '94:A', 95: '109:A'},
                {(69, 80), (69, 95), (80, 95)}),
]

@pytest.mark.parametrize("motif_possibilities,probe_coordinators, \
                            name_for_res,res_for_column,solution", testdata)

def test_check_actual_motif(motif_possibilities, probe_coordinators, 
                                name_for_res, res_for_column, solution):
    
    result = motif._check_actual_motif(motif_possibilities, 
                        probe_coordinators, name_for_res, res_for_column)
    
    assert solution == result