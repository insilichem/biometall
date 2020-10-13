import pytest
from pytest import approx
import numpy as np
from biometall.modules import grid

# Files: ["4dc8.pdb", "1hnn.pdb", "3rzy.pdb"]
testdata = [
    ([-4.3311, -6.3593, 15.3793], 29.6255, 107392, [-32.9523, -13.8912, 14.8772], [24.2901, 1.1726, 15.8814]),
    ([12.8607, 37.3031, 18.2772], 47.2562, 434407, [-34.3955, 37.3031, 18.2772], [60.1168, 37.3031, 18.2772]),
    ([8.1566, 10.4758, 13.9826], 38.4276, 229489, [-30.2710, 10.4758, 13.9826], [46.5842, 10.4758, 13.9826])
]
@pytest.mark.parametrize("centroids,radii,len_grid,first_point,last_point", testdata)

def test_grid(centroids, radii, len_grid, first_point, last_point):
    centroid = np.array(centroids)
    points = grid._grid(centroid, radii, 1.0)
    assert len(points) == len_grid
    assert first_point == approx(points[0], abs=1e-3)
    assert last_point == approx(points[-1], abs=1e-3)
