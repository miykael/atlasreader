import numpy as np
from mni_atlas_reader import atlas_reader


def test_get_vox_coord():
    affine = np.eye(4)
    coords = (10, 10, 10)
    assert atlas_reader.get_vox_coord(affine, coords) == [10, 10, 10]
