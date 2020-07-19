import pytest

from biometall.modules import motif


testdata = [
            ([('HIS', 'ASP'), ('HIS', 'GLU')],
                {'HIS', 'ASP', 'GLU'}, 
                {16: ['ASP', 'GLU'], 17: ['HIS','SER'], 18: ['TYR','HIS']}, 
                {16: {'ASP', 'GLU'}, 17: {'HIS'}, 18: {'HIS'}}),
            ([('CYS', 'SER'), ('CYS', 'MET')],
                {'CYS', 'SER', 'MET'}, 
                {16: ['ASP', 'GLU'], 17: ['HIS','SER'], 18: ['TYR','HIS']}, 
                {}),
            ([('HIS',)],
                {'HIS'}, 
                {16: ['ASP', 'GLU'], 17: ['HIS','SER'], 18: ['TYR','HIS']}, 
                {17: {'HIS'}, 18: {'HIS'}}),
]

@pytest.mark.parametrize("motif_possibilities,residues, \
                            probe_coordinators,solution", testdata)

def test_check_possible_mutations(motif_possibilities, residues, 
                                probe_coordinators, solution):
    
    result = motif._check_possible_mutations(motif_possibilities, residues, 
                                                probe_coordinators)
    
    assert solution == result