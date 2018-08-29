import os
import numpy as np
from atlasreader import atlasreader
from nilearn.datasets import fetch_neurovault_motor_task
import pytest

STAT_IMG = fetch_neurovault_motor_task().images[0]
EXAMPLE_COORDS = dict(
    affine=np.array([[1, 0, 0, -90],
                     [0, 1, 0, -150],
                     [0, 0, 1, -80],
                     [0, 0, 0, 1]]),
    coords=[
        dict(
            ijk=[0, 0, 0],
            xyz=[-90, -150, -80]
        ),
        dict(
            ijk=[[10, 10, 10], [100, 50, 100]],
            xyz=[[-80, -140, -70], [10, -100, 20]]
        ),
        dict(
            ijk=[[54, 32, 20], [82, 205, 38], [32, 51, 82]],
            xyz=[[-36, -118, -60], [-8, 55, -42], [-58, -99, 2]]
        )
    ]
)


def test_coords_transform():
    aff = EXAMPLE_COORDS['affine']
    for coords in EXAMPLE_COORDS['coords']:
        ijk, xyz = coords['ijk'], coords['xyz']
        assert atlasreader.coord_xyz_to_ijk(aff, xyz) == ijk
        assert atlasreader.coord_ijk_to_xyz(aff, ijk) == xyz
    with pytest.raises(ValueError):
        atlasreader.coord_xyz_to_ijk(aff, [[10, 10], [20, 30]])
    with pytest.raises(ValueError):
        atlasreader.coord_ijk_to_xyz(aff, [[10, 10], [20, 30]])


def test_get_statmap_info():
    # general integration test to check that min_distance works
    # this will take a little while since it's running it twice
    for min_distance in [None, 20]:
        cdf, pdf = atlasreader.get_statmap_info(STAT_IMG,
                                                atlas=['Harvard_Oxford',
                                                       'AAL'],
                                                min_distance=min_distance)


def test_create_output(tmpdir):
    # create mock data
    stat_img_name = os.path.basename(STAT_IMG)[:-7]

    # temporary output
    output_dir = tmpdir.mkdir('mni_test')
    atlasreader.create_output(STAT_IMG,
                              atlas=['Harvard_Oxford'],
                              outdir=output_dir)

    # test if output exists and if the key .csv and .png files were created
    assert os.path.exists(output_dir)
    assert len(os.listdir(output_dir)) > 0
    assert os.path.isfile(
        os.path.join(output_dir, '{}_clusters.csv'.format(stat_img_name))
    )
    assert os.path.isfile(
        os.path.join(output_dir, '{}_peaks.csv'.format(stat_img_name))
    )
    assert os.path.isfile(
        os.path.join(output_dir, '{}.png'.format(stat_img_name))
    )
