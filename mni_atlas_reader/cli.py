"""
Functions for command line interface to generate cluster / peak summary
"""
import argparse
import os.path as op
from mni_atlas_reader.atlas_reader import (_ATLASES, check_atlases,
                                           create_output)

_ACCEPTED_ATLASES = _ATLASES + [a.lower() for a in _ATLASES]


def _check_limit(num, limits=[0, 100]):
    """
    Ensures that `num` is within provided `limits`

    Parameters
    ----------
    num : float
        Number to assess
    limits : list, optional
        Lower and upper bounds that `num` must be between to be considered
        valid
    """
    lo, hi = limits
    num = float(num)

    if num < lo or num > hi:
        raise ValueError('Provided value {} is outside expected limits {}.'
                         .format(num, limits))

    return num


def _get_parser():
    """ Reads command line arguments and returns input specifications """
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=op.abspath, metavar='file',
                        help='The full or relative path to the statistical map'
                             'from which cluster information should be '
                             'extracted.')
    parser.add_argument('-a', '--atlas', type=str, default='all', nargs='+',
                        choices=_ACCEPTED_ATLASES + ['all'], metavar='atlas',
                        help='Atlas(es) to use for examining anatomical '
                             'delineation of clusters in provided statistical '
                             'map. Default: all available atlases.')
    parser.add_argument('-t', '--threshold', type=float, default=2,
                        dest='voxel_thresh', metavar='threshold',
                        help='Value threshold that voxels in provided file '
                             'must surpass in order to be considered in '
                             'cluster extraction.')
    parser.add_argument('-c', '--cluster', type=float, default=5,
                        dest='cluster_extent', metavar='extent',
                        help='Required number of contiguous voxels for a '
                             'cluster to be retained for analysis.')
    parser.add_argument('-p', '--probability', type=_check_limit, default=5,
                        dest='prob_thresh', metavar='threshold',
                        help='Threshold to consider when using a '
                             'probabilistic atlas for extracting anatomical '
                             'cluster locations. Value will apply to all '
                             'request probabilistic atlases, and should range '
                             'between 0 and 100.')
    parser.add_argument('-o', '--outdir', type=str, default=None,
                        dest='outdir', metavar='outdir',
                        help='Output directory for created files. If it is '
                             'not specified, then output files are created in '
                             'the same directory as the statistical map that '
                             'is provided.')

    return parser.parse_args()


def main():
    """
    The primary entrypoint for calling atlas reader via the command line

    All parameters are read via argparse, so this should only be called from
    the command line!
    """

    opts = _get_parser()
    create_output(opts.filename,
                  atlas=check_atlases(opts.atlas),
                  voxel_thresh=opts.voxel_thresh,
                  cluster_extent=opts.cluster_extent,
                  prob_thresh=opts.prob_thresh,
                  outdir=opts.outdir)


if __name__ == '__main__':
    raise RuntimeError('`mni_atlas_reader/cli.py` should not be run '
                       'directly. Please `pip install` mni_atlas_reader and '
                       'use the `mni_atlas_reader` command, instead.')
