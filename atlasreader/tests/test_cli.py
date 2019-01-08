import os
from nilearn.datasets import fetch_neurovault_motor_task
import subprocess

STAT_IMG = fetch_neurovault_motor_task().images[0]


def test_cli(tmpdir):
    stat_img_name = os.path.basename(STAT_IMG)[:-7]

    # temporary output
    output_dir = tmpdir.mkdir('mni_test')
    ret = subprocess.run(['atlasreader',
                          '--atlas', 'harvard_oxford', 'aal',
                          '--threshold', '1.96',
                          '--probability', '5',
                          '--mindist', '20',
                          '--outdir', output_dir,
                          STAT_IMG, '20'])
    assert ret.returncode == 0

    # test if output exists and if the key .csv and .png files were created
    assert output_dir.exists()
    assert len(output_dir.listdir()) > 0
    assert output_dir.join('{}_clusters.csv'.format(stat_img_name)).isfile()
    assert output_dir.join('{}_peaks.csv'.format(stat_img_name)).isfile()
    assert output_dir.join('{}.png'.format(stat_img_name)).isfile()
