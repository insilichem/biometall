import pytest
from pytest import approx
import numpy as np
from biometall.modules import grid

testdata = [
    ([[0,0,0], [0,0,1], [0,1,0], [0,1,1], [1,0,0], [1,0,1], [1,1,0], [1,1,1]], [0,0,0], 1.7321),
    ([[0,0,0], [0,0,1], [0,1,0], [0,1,1], [1,0,0], [1,0,1], [1,1,0], [1,1,1], [0.5,0.5,0.5]], [0.5,0.5,0.5], 0.8660),
    ([[0,0,0], [0,0,1], [0,1,0], [0,1,1], [1,0,0], [0.5,0.5,0.5], [1,0,1], [1,1,0], [1,1,1]], [0.5,0.5,0.5], 0.8660),
    ([[0.5,0.5,0.5], [0,0,0], [0,0,1], [0,1,0], [0,1,1], [1,0,0], [1,0,1], [1,1,0], [1,1,1]], [0.5,0.5,0.5], 0.8660),
    ([[0.0,0.0,0.0], [0.0,0.0,1.0], [0.0,0.0,0.5]], [0.0,0.0,0.5], 0.5),
    ([[-0.17,0.56,-0.11], [0.03,0.36,-0.31], [0.03,0.36,-0.11], [0.03,0.36, 0.09], [0.03,0.56,-0.31],
      [0.03,0.56,-0.11], [0.03,0.56, 0.09], [0.03,0.76,-0.31], [0.03,0.76,-0.11], [0.03,0.76, 0.09],
 	  [0.23,0.16,-0.11], [0.23,0.36,-0.31], [0.23,0.36,-0.11], [0.23,0.36, 0.09], [0.23,0.56,-0.51],
 	  [0.23,0.56,-0.31], [0.23,0.56,-0.11], [0.23,0.56, 0.09], [0.23,0.56, 0.29], [0.23,0.76,-0.31],
 	  [0.23,0.76,-0.11], [0.23,0.76, 0.09], [0.23,0.96,-0.11], [0.43,0.36,-0.31], [0.43,0.36,-0.11],
 	  [0.43,0.36, 0.09], [0.43,0.56,-0.31], [0.43,0.56,-0.11], [0.43,0.56, 0.09], [0.43,0.76,-0.31],
 	  [0.43,0.76,-0.11], [0.43,0.76, 0.09], [0.63,0.56,-0.11]], [0.23,0.56,-0.11], 0.4)
]
@pytest.mark.parametrize("probes,center,radius", testdata)

def test_calculate_center_and_radius(probes, center, radius):
    center_result, radius_result = grid._calculate_center_and_radius(probes)
    np_center = np.array(center)
    assert np_center == approx(center_result, abs=1e-3)
    assert radius_result == approx(radius, abs=1e-3)
