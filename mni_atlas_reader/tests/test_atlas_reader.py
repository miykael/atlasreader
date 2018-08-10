import os
import numpy as np
from mni_atlas_reader import atlas_reader
from nilearn.datasets import fetch_neurovault_motor_task

def test_get_vox_coord():
    affine = np.eye(4)
    coords = (10, 10, 10)
    assert atlas_reader.get_vox_coord(affine, coords) == [10, 10, 10]


def test_create_output(tmpdir):

    # create mock data
    motor_images = fetch_neurovault_motor_task()
    stat_img = motor_images.images[0]
    stat_img_name = os.path.basename(stat_img)[:-7]

    # temporary output
    output_dir = tmpdir.mkdir('mni_test')


    atlas_reader.create_output(stat_img, atlas=['Harvard_Oxford'],
                               outDir=output_dir)

    # test if output exists and if the key .csv and .png files were created
    assert os.path.exists(output_dir)
    assert len(os.listdir(output_dir)) > 0
    assert os.path.isfile(
        os.path.join(output_dir, '{}.csv'.format(stat_img_name))
    )
    assert os.path.isfile(
        os.path.join(output_dir, '{}.png'.format(stat_img_name))
    )
