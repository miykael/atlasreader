import os
import numpy as np
from atlasreader import atlasreader
import nibabel as nb
from nilearn.datasets import fetch_neurovault_motor_task
import pytest
import pandas as pd

STAT_IMG = fetch_neurovault_motor_task().images[0]
EXAMPLE_COORDS = dict(
    affine=np.array([[1, 0, 0, -90],
                     [0, 1, 0, -150],
                     [0, 0, 1, -80],
                     [0, 0, 0, 1]]),
    coords=[
        dict(
            ijk=[0, 0, 0],
            xyz=np.array([[-90, -150, -80]])
        ),
        dict(
            ijk=[[10, 10, 10], [100, 50, 100]],
            xyz=np.array([[-80, -140, -70], [10, -100, 20]])
        ),
        dict(
            ijk=[[54, 32, 20], [82, 205, 38], [32, 51, 82]],
            xyz=np.array([[-36, -118, -60], [-8, 55, -42], [-58, -99, 2]])
        )],
    bad_coords=[
        dict(
            ijk_in=np.array([80, 80, 80]),
            ijk_out=np.array([80, 80, 80])
        ),
        dict(
            ijk_in=np.array([-1, 0, 0]),
            ijk_out=np.array([0, 0, 0])
        ),
        dict(
            ijk_in=np.array([80, 80, 100]),
            ijk_out=np.array([0, 0, 0])
        )
    ],
    bounding_shape=np.array([90, 90, 90])
)
EXPECTED_TABLES = dict(
    cluster=np.array([[42,  -25,  58,  6.66003, 36936],
                      [-36, -25,  55, -6.63604, 15012],
                      [45,  -19,  16,  5.76538,  7722],
                      [-15, -52, -26,  6.11673,  7101],
                      [18,  -55, -23, -5.90086,  5184],
                      [-36, -19,  19, -5.00687,   648]]),
    peak=np.array([[42,  -25,  58,  7.94135, 36936],
                   [-36, -25,  55, -7.94144, 15012],
                   [45,  -19,  16,  7.94135,  7722],
                   [-15, -52, -26,  7.94135,  7101],
                   [18,  -55, -23, -7.94144,  5184],
                   [-36, -19,  19, -6.21808,   648]])
)


def test_get_atlases():
    for atlas in atlasreader._ATLASES:
        a = atlasreader.get_atlas(atlas, cache=False)
        assert all(hasattr(a, k) for k in ['atlas', 'image', 'labels'])
    with pytest.raises(ValueError):
        atlasreader.get_atlas('not_an_atlas')


def test_check_atlases():
    atlases = atlasreader.check_atlases('all')
    assert len(atlases) == len(atlasreader._ATLASES)
    atlases = atlasreader.check_atlases('default')
    assert len(atlases) == len(atlasreader._DEFAULT)
    atlases = atlasreader.check_atlases(['aal', 'destrieux'])
    assert atlasreader.check_atlases(atlases) == atlases
    assert atlasreader.check_atlases(atlases[0]) == atlases[0]


def test_coords_transform():
    aff = EXAMPLE_COORDS['affine']
    for coords in EXAMPLE_COORDS['coords']:
        ijk, xyz = coords['ijk'], coords['xyz']
        assert np.all(atlasreader.coord_xyz_to_ijk(aff, xyz) == ijk)
        assert np.all(atlasreader.coord_ijk_to_xyz(aff, ijk) == xyz)
    with pytest.raises(ValueError):
        atlasreader.coord_xyz_to_ijk(aff, [[10, 10], [20, 30]])
    with pytest.raises(ValueError):
        atlasreader.coord_ijk_to_xyz(aff, [[10, 10], [20, 30]])


def test_bounding_box_check():
    for coords in EXAMPLE_COORDS['bad_coords']:
        ijk_out = atlasreader.check_atlas_bounding_box(
            coords['ijk_in'], EXAMPLE_COORDS['bounding_shape'])
        assert np.all(ijk_out == coords['ijk_out'])


def test_get_statmap_info():
    # general integration test to check that min_distance works
    # this will take a little while since it's running it twice
    stat_img = nb.load(STAT_IMG)
    for min_distance in [None, 20]:
        cdf, pdf = atlasreader.get_statmap_info(stat_img,
                                                cluster_extent=20,
                                                atlas=['Harvard_Oxford',
                                                       'AAL'],
                                                min_distance=min_distance)

    # test that empty image return empty dataframes
    zero_img = nb.Nifti1Image(np.zeros(stat_img.shape), stat_img.affine,
                              header=stat_img.header)
    cdf, pdf = atlasreader.get_statmap_info(zero_img,
                                            cluster_extent=20)
    assert len(cdf) == 0
    assert len(pdf) == 0


def test_read_atlas_preaks():
    # Load a correct atlas
    atlasreader.read_atlas_peak('aal', [10, 10, 10])
    # Load a list of atlases
    with pytest.raises(ValueError):
        atlasreader.read_atlas_peak(2*['aal'], [10, 10, 10])
    # Load 'all' atlas
    with pytest.raises(ValueError):
        atlasreader.read_atlas_peak('all', [10, 10, 10])


def test_process_image():
    stat_img = nb.load(STAT_IMG)
    # check that defaults for processing image work
    img = atlasreader.process_img(stat_img, cluster_extent=20)
    assert isinstance(img, nb.Nifti1Image)
    # check that one-sided thresholding works
    img = atlasreader.process_img(stat_img, direction='neg', cluster_extent=20)
    assert isinstance(img, nb.Nifti1Image)
    # check that negative voxel threshold works
    img = atlasreader.process_img(stat_img, cluster_extent=20,
                                  voxel_thresh=-10)
    assert isinstance(img, nb.Nifti1Image)
    # check that setting cluster extent too high still returns an image
    img = atlasreader.process_img(stat_img, cluster_extent=5000)
    assert np.allclose(img.get_data(), 0)
    # ensure empty image --> empty image
    zero_img = nb.Nifti1Image(np.zeros(stat_img.shape), stat_img.affine,
                              header=stat_img.header)
    img = atlasreader.process_img(zero_img, cluster_extent=20)
    assert img.shape == zero_img.shape + (1,)
    assert np.allclose(img.get_data(), 0)


def test_create_output(tmpdir):
    # create mock data
    stat_img_name = os.path.basename(STAT_IMG)[:-7]

    # temporary output
    output_dir = tmpdir.mkdir('mni_test')
    atlasreader.create_output(STAT_IMG, cluster_extent=20,
                              voxel_thresh=7,
                              atlas=['Harvard_Oxford'],
                              outdir=output_dir)

    # test if output exists and if the key .csv and .png files were created
    assert output_dir.exists()
    assert len(output_dir.listdir()) > 0
    assert output_dir.join('{}_clusters.csv'.format(stat_img_name)).isfile()
    assert output_dir.join('{}_peaks.csv'.format(stat_img_name)).isfile()
    assert output_dir.join('{}.png'.format(stat_img_name)).isfile()


def test_plotting(tmpdir):
    """Test functionality of kwarg implementation"""

    # temporary output
    output_dir = tmpdir.mkdir('mni_test')

    # overwrite some default params
    atlasreader.create_output(STAT_IMG, cluster_extent=20,
                              voxel_thresh=7,
                              atlas=['Harvard_Oxford'],
                              outdir=output_dir,
                              glass_plot_kws={'display_mode': 'ortho'},
                              stat_plot_kws={'black_bg': False})

    # add new parameter not already set by default
    atlasreader.create_output(STAT_IMG, cluster_extent=20,
                              voxel_thresh=7,
                              atlas=['Harvard_Oxford'],
                              outdir=output_dir,
                              glass_plot_kws={'alpha': .4})


def test_table_output(tmpdir):
    # create mock data
    stat_img_name = os.path.basename(STAT_IMG)[:-7]

    # temporary output
    output_dir = tmpdir.mkdir('mni_test')
    atlasreader.create_output(STAT_IMG, cluster_extent=20,
                              voxel_thresh=4, atlas='default',
                              outdir=output_dir)

    # test if output tables contain expected output
    df = pd.read_csv(output_dir.join('{}_clusters.csv'.format(stat_img_name)))
    assert np.allclose(df[df.keys()[1:6]].values, EXPECTED_TABLES['cluster'])
    df = pd.read_csv(output_dir.join('{}_peaks.csv'.format(stat_img_name)))
    assert np.allclose(df[df.keys()[1:6]].values, EXPECTED_TABLES['peak'])
