import pytest
from biometall.modules import grid

testdata = [
    (19870,437,1,268435456,9936),
    (19870,437,4,268435456,2484),
    (19870,437,8,268435456,1325),
    (19870,437,1,2147483648,19871),
    (19870,437,4,2147483648,4968),
    (19870,437,8,2147483648,2484),
    (19870,437,1,4294967296,19871),
    (19870,437,4,4294967296,4968),
    (19870,437,8,4294967296,2484),
]
@pytest.mark.parametrize("len_grid,len_protein,n_cores,optimize_to,chunk_size", testdata)

def test_chunk_size(len_grid, len_protein, n_cores, optimize_to, chunk_size):
    chunk_size_result = grid._chunk_size(len_grid, len_protein, n_cores, optimize_to)
    assert chunk_size == chunk_size_result
