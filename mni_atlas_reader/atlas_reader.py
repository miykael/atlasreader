"""
Primary functions of mni_atlas_reader
"""
import argparse
import os
import os.path as op
from pkg_resources import resource_filename
import nibabel as nb
from nilearn.plotting import plot_glass_brain, plot_stat_map
import numpy as np
import pandas as pd
from scipy.ndimage import label
from sklearn.utils import Bunch


_ATLASES = [
    'AAL', 'Desikan-Killiany', 'Destrieux', 'Harvard_Oxford', 'Juelich',
    'Neuromorphometrics',
]
_ACCEPTED_ATLASES = _ATLASES + [a.lower() for a in _ATLASES]


def get_atlas(atlastype):
    """
    Gets `atlastype` image and corresponding label file from package resources

    Parameters
    ----------
    atlastype : str
        Name of atlas to query

    Returns
    -------
    info : sklearn.utils.Bunch
        image : Niimg_like
            ROI image loaded with integer-based labels indicating parcels
        labels : pandas.core.data.DataFrame
            Dataframe with columns ['index', 'name'] matching region IDs in
            `image` to anatomical `name`
    """
    # we accept uppercase atlases via argparse but our filenames are lowercase
    atlastype = atlastype.lower()

    # get the path to atlas + label files shipped with package
    # resource_filename ensures that we're getting the correct path
    data_dir = resource_filename('mni_atlas_reader', 'data/atlases')
    atlas_path = op.join(data_dir, 'atlas_{0}.nii.gz'.format(atlastype))
    label_path = op.join(data_dir, 'labels_{0}.csv'.format(atlastype))

    # return loaded filenames (we should only have to call this once!)
    return Bunch(image=nb.load(atlas_path), labels=pd.read_csv(label_path))


def _check_coord_inputs(coords):
    """
    Confirms `coords` are appropriate shape for coordinate transform

    Parameters
    ----------
    coords : array-like

    Returns
    -------
    coords : (4 x N) numpy.ndarray
    """
    coords = np.atleast_2d(coords).T
    if 3 not in coords.shape:
        raise ValueError('Input coordinates must be of shape (3 x N). '
                         'Provided coordinate shape: {}'.format(coords.shape))
    if coords.shape[0] != 3:
        coords = coords.T
    # add constant term to coords to make 4 x N
    coords = np.row_stack([coords, np.ones_like(coords[0])])
    return coords


def coord_ijk_to_xyz(affine, coords):
    """
    Converts voxel `coords` in cartesian space to `affine` space

    Parameters
    ----------
    affine : (4, 4) array-like
        Affine matrix
    coords : (N,) list of list
        Image coordinate values, where each entry is a length three list of int
        denoting ijk coordinates in cartesian space

    Returns
    ------
    xyz : (N,) list of list
        Provided `coords` in `affine` space, where each entry is a length three
        list of float denoting xyz coordinates
    """
    coords = _check_coord_inputs(coords)
    mni_coords = np.dot(affine, coords)[:3].T
    return mni_coords.squeeze().tolist()


def coord_xyz_to_ijk(affine, coords):
    """
    Converts voxel `coords` in `affine` space to cartesian space

    Parameters
    ----------
    affine : (4, 4) array-like
        Affine matrix
    coords : (N,) list of list
        Image coordinate values, where each entry is a length three list of int
        denoting xyz coordinates in `affine` space

    Returns
    ------
    ijk : (N,) list of list
        Provided `coords` in cartesian space, where each entry is a length
        three list of float denoting ijk coordinates
    """
    coords = _check_coord_inputs(coords)
    vox_coords = np.linalg.solve(affine, coords)[:3].T.astype(int)
    return vox_coords.squeeze().tolist()


def get_label(atlastype, label_id):
    """
    Gets anatomical name of `label_id` in `atlastype`

    Parameters
    ----------
    atlastype : str
        Name of atlas to use
    label_id : int
        Numerical ID representing label

    Returns
    ------
    label : str
        Neuroanatomical region of `label_id` in `atlastype`
    """
    labels = get_atlas(atlastype).labels
    try:
        return labels.query('index == {}'.format(label_id)).name.iloc[0]
    except IndexError:
        return 'no_label'


def get_clusters(data, min_extent=5):
    """
    Extracts clusters from statistical map `data`

    Parameters
    ----------
    data : (X, Y, Z) array-like
        Thresholded data from statistical map
    min_extent : int, optional
        Minimum number of voxels required to consider a cluster. Default: 5

    Returns
    ------
    clusters : numpy.ndarray
        Array of numerically labelled clusters
    nclusters : int
        Number of clusters in `data`
    """
    clusters, nclusters = label(data)
    for idx in range(1, nclusters + 1):
        if np.sum(clusters == idx) < min_extent:
            clusters[clusters == idx] = 0
    nclusters = len(np.setdiff1d(np.unique(clusters), [0]))
    return clusters, nclusters


def get_peak_coords(img, affine, data):
    """
    Gets MNI coordinates of peak voxels within each cluster of `data`

    Parameters
    ----------
    img : (X, Y, Z) array-like
        Array of numerically labelled clusters
    affine : (4, 4) array-like
        Affine matrix
    data : (X, Y, Z) array-like
        Thresholded data from statistical map

    Returns
    ------
    coordinates : list of lists
        Coordinates of peak voxels in `data`
    """
    coords = []
    clusters = np.setdiff1d(np.unique(img.ravel()), [0])
    clust_size = []
    maxcoords = []
    for lab in clusters:
        clust_size.append(np.sum(img == lab))  # get cluster size
        maxval = np.max(data[img == lab])      # find peak value of cluster
        maxidx = np.nonzero(np.multiply(data, img == lab) == maxval)
        maxcoords.append([m[0] for m in maxidx])

    maxcoords = np.asarray(maxcoords)
    maxcoords = maxcoords[np.argsort(clust_size)[::-1], :]
    for i, lab in enumerate(clusters[np.argsort(clust_size)[::-1]]):
        coords.append(coord_ijk_to_xyz(affine, maxcoords[i]))
    return coords


def get_cluster_coords(cluster, affine):
    """
    Get `affine` coordinates of voxels in `cluster`

    Parameters
    ----------
    cluster : (X, Y, Z) array-like
        Boolean mask of a single cluster
    affine : (4, 4) array-like
        Affine matrix

    Returns
    ------
    coordinates : list of numpy.ndarray
        List of coordinates for voxels in cluster
    """
    coords_vox = np.rollaxis(np.array(np.where(cluster)), 1)
    coords = [coord_ijk_to_xyz(affine, c) for c in coords_vox]
    return coords


def read_atlas_peak(atlastype, coordinate, prob_thresh=5):
    """
    Returns label of `coordinate` from corresponding `atlastype`

    If `atlastype` is probabilistic, `prob_thresh` determines (in percentage
    units) the threshold to apply before getting label of `coordinate`

    Parameters
    ----------
    atlastype : str
        Name of atlas to use
    coordinate : list of float
        x, y, z MNI coordinates of voxel
    prob_thresh : [0, 100] int, optional
        Probability (percentage) threshold to apply if `atlastype` is
        probabilistic

    Returns
    ------
    label : str or list of lists
        If `atlastype` is deterministic, this is the corresponding atlas label
        of `coordinate`. If `atlastype` is probabilistic, this is a list of
        lists where each entry denotes the probability and corresponding label
        of `coordinate`.
    """

    atlas = get_atlas(atlastype).image

    # get atlas data and affine matrix
    data = atlas.get_data()
    affine = atlas.affine

    # get voxel index
    voxID = coord_xyz_to_ijk(affine, coordinate)

    # get label information
    if atlastype.lower() in ['juelich', 'harvard_oxford']:
        probs = data[voxID[0], voxID[1], voxID[2]]
        probs[probs < prob_thresh] = 0
        idx = np.where(probs)[0]

        # sort list by probability
        idx = idx[np.argsort(probs[idx])][::-1]

        # get probability and label names
        probLabel = []
        for i in idx:
            label = get_label(atlastype, i)
            probLabel.append([probs[i], label])

        # if no labels found
        if probLabel == []:
            probLabel = [[0, 'no_label']]

        return probLabel

    else:
        labelID = int(data[voxID[0], voxID[1], voxID[2]])
        label = get_label(atlastype, labelID)
        return label


def read_atlas_cluster(atlastype, cluster, affine, prob_thresh=5):
    """
    Returns label of `cluster` from corresponding `atlastype`

    If `atlastype` is probabilistic, `prob_thresh` determines (in percentage
    units) the threshold to apply before getting label of `coordinate`

    Parameters
    ----------
    atlastype : str
        Name of atlas to use
    cluster : (X, Y, Z) array-like
        Boolean mask of cluster
    coordinate : list of float
        x, y, z MNI coordinates of voxel
    prob_thresh : [0, 100] int, optional
        Probability (percentage) threshold to apply if `atlastype` is
        probabilistic

    Returns
    ------
    segments : list of lists
        Where each entry is of the form [probability, label] denoting the
        extent to which `cluster` overlaps with region `label` in `atlastype`
    """
    atlas = get_atlas(atlastype).image

    # get atlas data and affine matrix
    atlas_data = atlas.get_data()
    atlas_affine = atlas.affine

    # get coordinates of each voxel in cluster
    coords = get_cluster_coords(cluster, affine)

    # get voxel indexes
    voxIDs = [coord_xyz_to_ijk(atlas_affine, c) for c in coords]

    # get label information
    if atlastype.lower() in ['juelich', 'harvard_oxford']:
        labelIDs = [np.argmax(atlas_data[v[0], v[1], v[2]]) if np.sum(
            atlas_data[v[0], v[1], v[2]]) != 0 else -1 for v in voxIDs]
    else:
        labelIDs = [int(atlas_data[v[0], v[1], v[2]]) for v in voxIDs]

    unique_labels = np.unique(labelIDs)
    labels = np.array([get_label(atlastype, u) for u in unique_labels])
    N = float(len(labelIDs))
    percentage = np.array(
        [100 * np.sum(labelIDs == u) / N for u in unique_labels])

    sortID = np.argsort(percentage)[::-1]

    return [[percentage[s], labels[s]] for s in sortID if
            percentage[s] >= prob_thresh]


def get_peak_info(coord, atlastype='all', prob_thresh=5):
    """
    Gets region and probability information for `coord` in `atlastype`

    Parameters
    ----------
    coordinate : list of float
        x, y, z MNI coordinates of voxel
    atlastype : str or list, optional
        Name of atlas(es) to use. Default: 'all'
    prob_thresh : [0, 100] int, optional
        Probability (percentage) threshold to apply to `atlastype` if it is
        probabilistic

    Returns
    -------
    peakinfo : list of lists
        Where each entry contains the atlas name and probability information
        of region labels associated with `coordinate` in the atlas
    """
    peakinfo = []
    if isinstance(atlastype, str):
        atlastype = [atlastype]

    if 'all' not in atlastype:
        for atypes in atlastype:
            segment = read_atlas_peak(atypes, coord, prob_thresh)
            peakinfo.append([atypes, segment])
    else:
        for atypes in _ATLASES:
            segment = read_atlas_peak(atypes, coord, prob_thresh)
            peakinfo.append([atypes, segment])

    return peakinfo


def get_cluster_info(cluster, affine, atlastype='all', prob_thresh=5):
    """
    Gets region and probability information for `cluster` in `atlastype`

    Parameters
    ----------
    cluster : (X, Y, Z) array-like
        Boolean mask of cluster
    affine : (4, 4) array-like
        Affine matrix mapping `cluster` to MNI space
    atlastype : str or list, optional
        Name of atlas(es) to use. Default: 'all'
    prob_thresh : [0, 100] int, optional
        Probability (percentage) threshold to apply to `atlastype`, if it is
        probabilistic

    Returns
    ------
    clusterinfo : list of lists
        Where each entry contains the atlas name and probability information of
        region labels associated with `cluster` in the atlas
    """

    clusterinfo = []

    if isinstance(atlastype, str):
        atlastype = [atlastype]

    if 'all' not in atlastype:
        for atypes in atlastype:
            segment = read_atlas_cluster(atypes, cluster, affine, prob_thresh)
            clusterinfo.append([atypes, segment])
    else:
        for atypes in _ATLASES:
            segment = read_atlas_cluster(atypes, cluster, affine, prob_thresh)
            clusterinfo.append([atypes, segment])

    return clusterinfo


def create_output(filename, atlas='all', voxel_thresh=1.96,
                  cluster_extent=20, prob_thresh=5, outdir=None):
    """
    Performs full cluster analysis on `filename`

    Generates output table containing information on each cluster in `filename`
    including: (1) number of voxels in cluster, (2) average activation across
    voxels, (3) MNI coordinates of peak voxel, and (4) neuroanatomical location
    of peak voxel based on specified `atlas`.

    In addition, screenshots of statistical maps are separately created for
    each cluster in which cross hairs are focused on the cluster peak voxel.

    Parameters
    ----------
    filename : str
        Path to input statistical map
    atlas : str, optional
        Name of atlas(es) to consider for cluster analysis. Default: 'all'
    voxel_thresh : int, optional
        Threshold applied to `filename` before performing cluster analysis.
        Default: 1.96
    cluster_extent : int, optional
        Minimum number of contiguous voxels required to consider a cluster in
        `filename`. Default: 20
    prob_thresh : int, optional
        Probability (percentage) threshold to apply to `atlas`, if it is
        probabilistic. Default: 5
    outdir : str or None, optional
        Path to desired output directory. If None, generated files will be
        saved to the same folder as `filename`. Default: None
    """

    fname = os.path.abspath(filename)

    # set up output directory
    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)
        savedir = outdir
    else:
        savedir = os.path.dirname(fname)

    # set up output filename, which is the same name as input w/o extension
    out_fname = os.path.basename(fname).split('.')[0]

    # get data from NIfTI file
    img = nb.load(filename)
    imgdata = img.get_data()
    if len(imgdata.shape) != 3:
        imgdata = imgdata[:, :, :, 0]

    # get top x-% of voxels if voxel_thresh is negative
    if voxel_thresh < 0:
        voxel_thresh = np.percentile(
            np.abs(imgdata[imgdata != 0]), (100 + voxel_thresh))

    # get clusters from data
    clusters, nclusters = get_clusters(
        np.abs(imgdata) > voxel_thresh, min_extent=cluster_extent)

    # clean img data
    imgdata[clusters == 0] = 0
    new_image = nb.Nifti1Image(imgdata, img.affine, img.header)

    # plot glass brain
    color_max = np.array([imgdata.min(), imgdata.max()])
    color_max = np.abs(np.min(color_max[color_max != 0]))
    glass_file = os.path.join(savedir, '{}.png'.format(out_fname))
    try:
        plot_glass_brain(new_image, vmax=color_max, threshold='auto',
                         display_mode='lyrz', black_bg=True,
                         plot_abs=False, colorbar=True,
                         output_file=glass_file)
    except ValueError:
        plot_glass_brain(new_image, vmax=color_max, threshold='auto',
                         black_bg=True, plot_abs=False, colorbar=True,
                         output_file=glass_file)

    # get coordinates of peaks
    coords = get_peak_coords(clusters, img.affine, np.abs(imgdata))

    # get peak and cluster information
    peak_summary = []
    peak_value = []
    cluster_summary = []
    cluster_mean = []
    volume_summary = []

    for n, coord in enumerate(coords):
        peakinfo = get_peak_info(
            coord, atlastype=atlas, prob_thresh=prob_thresh)
        peak_summary.append([p[1] if type(p[1]) != list else '; '.join(
            ['% '.join(e) for e in np.array(p[1])]) for p in peakinfo])
        voxID = coord_xyz_to_ijk(img.affine, coord)
        peak_value.append(imgdata[voxID[0], voxID[1], voxID[2]])

        idx = coord_xyz_to_ijk(img.affine, coord)
        clusterID = clusters[idx[0], idx[1], idx[2]]
        clusterinfo = get_cluster_info(
            clusters == clusterID,
            img.affine,
            atlastype=atlas,
            prob_thresh=prob_thresh)
        cluster_summary.append(['; '.join(
            ['% '.join([str(round(e[0], 2)), e[1]]) for e in coord[1]]) for
             coord in clusterinfo])
        cluster_mean.append(imgdata[clusters == clusterID].mean())

        voxel_volume = int(img.header['pixdim'][1:4].prod())
        volume_summary.append(np.sum(clusters == clusterID) * voxel_volume)

    # write output .csv file
    header = [p[0] for p in peakinfo]
    cluster_fname = os.path.join(savedir, '{}_clusters.csv'.format(out_fname))
    with open(cluster_fname, 'w') as f:
        f.writelines(','.join(
            ['ClusterID', 'Peak_Location', 'Cluster_Mean', 'Volume'] + header)
            + '\n'
        )

        for i, c in enumerate(cluster_summary):
            f.writelines(
                ','.join(['Cluster%.02d' % (i + 1), '_'.join(
                    [str(xyz) for xyz in coords[i]]), str(cluster_mean[i]),
                     str(volume_summary[i])] + c) + '\n')

        f.writelines('\n')

    peaks_fname = os.path.join(savedir, '{}_peaks.csv'.format(out_fname))
    with open(peaks_fname, 'w') as f:
        f.writelines(','.join(
            ['PeakID', 'Peak_Location', 'Peak_Value', 'Volume'] + header)
            + '\n'
        )

        for i, p in enumerate(peak_summary):
            f.writelines(
                ','.join(['Peak%.02d' % (i + 1), '_'.join(
                    [str(xyz) for xyz in coords[i]]), str(peak_value[i]),
                     str(volume_summary[i])] + p) + '\n')

        f.writelines('\n')

    # get template image for plotting cluster maps
    bgimg = nb.load(
        resource_filename(
            'mni_atlas_reader',
            'data/templates/MNI152_T1_1mm_brain.nii.gz'
        )
    )
    # plot clusters
    for idx, coord in enumerate(coords):
        cluster_name = '{}_cluster{:02d}'.format(out_fname, idx + 1)
        out_cluster_file = os.path.join(savedir, '{}.png'.format(cluster_name))

        try:
            plot_stat_map(new_image, vmax=color_max,
                          colorbar=True, title=cluster_name,
                          threshold=voxel_thresh, draw_cross=True,
                          black_bg=True, symmetric_cbar=True,
                          output_file=out_cluster_file,
                          bg_img=bgimg, cut_coords=coord, display_mode='ortho')
        except ValueError:
            plot_stat_map(new_image, vmax=color_max,
                          colorbar=True, title=cluster_name,
                          threshold=voxel_thresh, draw_cross=True,
                          black_bg=True, symmetric_cbar=True,
                          output_file=out_cluster_file)


def check_limit(num, limits=[0, 100]):
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
    parser.add_argument('-p', '--probability', type=check_limit, default=5,
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
                  atlas=opts.atlas,
                  voxel_thresh=opts.voxel_thresh,
                  cluster_extent=opts.cluster_extent,
                  prob_thresh=opts.prob_thresh,
                  outdir=opts.outdir)


if __name__ == '__main__':
    raise RuntimeError('`mni_atlas_reader/atlas_reader.py` should not be run '
                       'directly. Please `pip install` mni_atlas_reader and '
                       'use the `mni_atlas_reader` command, instead.')
