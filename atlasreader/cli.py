"""
Functions for command line interface to generate cluster / peak summary
"""
import argparse
import os.path as op
from atlasreader.atlasreader import (_ATLASES, check_atlases, create_output,
                                     read_atlas_peak)


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


def _atlasreader_parser():
    """ Reads command line arguments and returns input specifications """
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=op.abspath, metavar='file',
                        help='The full or relative path to the statistical '
                             'map from which cluster information should be '
                             'extracted.')
    parser.add_argument('cluster_extent', type=int, metavar='min_cluster_size',
                        help='Number of contiguous voxels required for a '
                             'cluster to be considered for analysis.')
    parser.add_argument('-a', '--atlas', type=str.lower, default='default',
                        nargs='+', choices=_ATLASES + ['all', 'default'],
                        metavar='atlas',
                        help='Atlas(es) to use for examining anatomical '
                             'delineation of clusters in provided statistical '
                             'map. Default: AAL, Desikan-Killiany & '
                             'Harvard-Oxford.')
    parser.add_argument('-t', '--threshold', type=float, default=2,
                        dest='voxel_thresh', metavar='threshold',
                        help='Value threshold that voxels in provided file '
                             'must surpass in order to be considered in '
                             'cluster extraction. Default: 2')
    parser.add_argument('-p', '--probability', type=_check_limit, default=5,
                        dest='prob_thresh', metavar='threshold',
                        help='Threshold to consider when using a '
                             'probabilistic atlas for extracting anatomical '
                             'cluster locations. Value will apply to all '
                             'request probabilistic atlases, and should range '
                             'between 0 and 100. Default: 5')
    parser.add_argument('-o', '--outdir', type=str, default=None,
                        dest='outdir', metavar='outdir',
                        help='Output directory for created files. If it is '
                             'not specified, then output files are created in '
                             'the same directory as the statistical map that '
                             'is provided. Default: None')
    parser.add_argument('-d', '--mindist', type=float, default=None,
                        dest='min_distance', metavar='distance',
                        help='If specified, the program will attempt to find '
                             'subpeaks within detected clusters, rather than '
                             'a single peak per cluster. The specified value '
                             'will determine the minimum distance (in mm) '
                             'required between subpeaks. Default: None')
    return parser.parse_args()


def atlasreader_main():
    """
    The primary entrypoint for calling atlas reader via the command line

    All parameters are read via argparse, so this should only be called from
    the command line!
    """

    opts = _atlasreader_parser()
    create_output(opts.filename,
                  atlas=check_atlases(opts.atlas),
                  voxel_thresh=opts.voxel_thresh,
                  cluster_extent=opts.cluster_extent,
                  prob_thresh=opts.prob_thresh,
                  outdir=opts.outdir,
                  min_distance=opts.min_distance)


def _queryatlas_parser():
    """ Reads command line arguments and returns input specifications """

    parser = argparse.ArgumentParser()
    parser.add_argument('coordinate', nargs=3, type=float,
                        help='The coordinate (in MNI space) for which to '
                             'query atlases. Should be provided as space-'
                             'delimited x y z coordinates, e.g., -10.5 -30.0 '
                             '15.3.')
    parser.add_argument('-a', '--atlas', type=str.lower, default='all',
                        nargs='+', choices=_ATLASES + ['all'], metavar='atlas',
                        help='Atlas(es) to use for examining anatomical '
                             'delineation of clusters in provided statistical '
                             'map. Default: all available atlases.')
    parser.add_argument('-p', '--probability', type=_check_limit, default=5,
                        dest='prob_thresh', metavar='threshold',
                        help='Threshold to consider when using a '
                             'probabilistic atlas for extracting anatomical '
                             'cluster locations. Value will apply to all '
                             'requested probabilistic atlases, and should '
                             'range between 0 and 100. Default: 5')
    return parser.parse_args()


def queryatlas_main():
    """
    The primary entrypoint for querying atlases via the command line

    All parameters are read via argparse, so this should only be called from
    the command line!
    """

    opts = _queryatlas_parser()
    print('{0:<25} {1:<25}\n{2:<25} {2:<25}'.format('Atlas', 'Label', '='*10))
    for atlas in check_atlases(opts.atlas):
        label = read_atlas_peak(atlas, opts.coordinate, opts.prob_thresh)
        if isinstance(label, list):
            label = '\n{}'.format(' '*26).join(
                ['{:>2}% {}'.format(*e) for e in label]
            )
        print('{:<25} {}'.format(atlas.atlas, label))


if __name__ == '__main__':
    raise RuntimeError('`atlasreader/cli.py` should not be run '
                       'directly. Please `pip install` atlasreader and '
                       'use the `atlasreader` or `queryatlas` commands, '
                       'instead.')
