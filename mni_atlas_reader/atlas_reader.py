"""
Primary functions of mni_atlas_reader
"""
import os
import os.path as op
from pkg_resources import resource_filename
import warnings
import nibabel as nb
from nilearn import image, plotting
from nilearn.regions import connected_regions
from nilearn._utils import check_niimg
import numpy as np
import pandas as pd
from scipy import ndimage
from sklearn.utils import Bunch


_ATLASES = [
    'AAL', 'Desikan_Killiany', 'Destrieux', 'Harvard_Oxford', 'Juelich',
    'Neuromorphometrics',
]


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


def clusterize_img(stat_img, cluster_extent=20):
    """
    Extracts clusters from statistical map `stat_img`

    Parameters
    ----------
    stat_img : Niimg_like object
        Thresholded statistical map image
    cluster_extent : int, optional
        Minimum number of voxels required to consider a cluster. Default: 20

    Returns
    ------
    cluster_img : Nifti1Image
        4D image of brain regions, where each volume is a distinct cluster
    """

    stat_img = check_niimg(stat_img)
    min_region_size = cluster_extent * np.prod(stat_img.header.get_zooms())
    cluster_img = connected_regions(stat_img,
                                    min_region_size=min_region_size,
                                    extract_type='connected_components')[0]

    return cluster_img


def get_peak_coords(cluster_img):
    """
    Gets MNI coordinates of peak voxels within each cluster of `cluster_img`

    Parameters
    ----------
    cluster_img : 4D-niimg_like
        4D image of brain regions, where each volume is a separated cluster

    Returns
    ------
    coordinates : list of lists
        Coordinates of peak voxels in `cluster_img`
    """

    # check cluster image and make it 4D, if not already
    cluster_img = check_niimg(cluster_img, atleast_4d=True)

    # create empty arrays to hold cluster size + peak coordinates
    clust_size = np.zeros(cluster_img.shape[-1])
    maxcoords = np.zeros((cluster_img.shape[-1], 3))

    # iterate through clusters and get info
    for n, cluster in enumerate(image.iter_img(cluster_img)):
        cluster = np.abs(cluster.get_data())
        clust_size[n] = np.sum(cluster != 0)
        maxcoords[n] = ndimage.center_of_mass(cluster == cluster.max())

    # sort peak coordinates by cluster size
    maxcoords = np.floor(maxcoords)[np.argsort(clust_size)[::-1]]

    # convert coordinates to MNI space
    coords = coord_ijk_to_xyz(cluster_img.affine, maxcoords)

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
    coords = coord_ijk_to_xyz(affine, coords_vox)
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
    # probabilistic atlas is requested
    if atlastype.lower() in ['juelich', 'harvard_oxford']:
        probs = data[voxID[0], voxID[1], voxID[2]]
        probs[probs < prob_thresh] = 0
        idx = np.where(probs)[0]

        # if no labels found
        if len(idx) == 0:
            return [[0, 'no_label']]

        # sort list by probability
        idx = idx[np.argsort(probs[idx])][::-1]

        # get probability and label names
        probLabel = [[probs[i], get_label(atlastype, i)] for i in idx]

        return probLabel
    # non-probabilistic atlas is requested
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
    affine : (4, 4) array-like
        Affine matrix
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
    voxIDs = coord_xyz_to_ijk(atlas_affine, coords)

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
    for atypes in check_atlases(atlastype):
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
    for atypes in check_atlases(atlastype):
        segment = read_atlas_cluster(atypes, cluster, affine, prob_thresh)
        clusterinfo.append([atypes, segment])

    return clusterinfo


def check_atlases(atlastype):
    """
    Converts `atlastype` to list and expands 'all', if present

    Parameters
    ----------
    atlastype : str or list
        Name of atlas(es) to use

    Returns
    -------
    atlases : list
        Names of atlas(es) to use
    """
    if isinstance(atlastype, str):
        atlastype = [atlastype]

    if 'all' not in atlastype:
        return atlastype
    else:
        return _ATLASES.copy()


def process_img(stat_img, voxel_thresh=1.96, cluster_extent=20):
    """
    Parameters
    ----------
    stat_img : Niimg_like object
        Thresholded statistical map image
    voxel_thresh : int, optional
        Threshold to apply to `stat_img`. If a negative number is provided a
        percentile threshold is used instead, where the percentile is
        determined by the equation `100 - voxel_thresh`. Default: 1.96
    cluster_extent : int, optional
        Minimum number of voxels required to consider a cluster. Default: 20

    Returns
    -------
    cluster_img : Nifti1Image
        4D image of brain regions, where each volume is a distinct cluster
    """
    # get input data image
    stat_img = image.index_img(check_niimg(stat_img, atleast_4d=True), 0)

    # threshold image
    if voxel_thresh < 0:
        voxel_thresh = '{}%'.format(100 + voxel_thresh)
    thresh_img = image.threshold_img(stat_img, threshold=voxel_thresh)

    # extract clusters
    cluster_img = clusterize_img(thresh_img, cluster_extent=cluster_extent)

    return cluster_img


def get_statmap_info(stat_img, atlas='all', voxel_thresh=1.96,
                     cluster_extent=20, prob_thresh=5):
    """
    Extract peaks and cluster information from `clust_img` for `atlas`

    Parameters
    ----------
    cluster_img : Niimg_like
        4D image of brain regions, where each volume is a distinct cluster
    atlas : str or list, optional
        Name of atlas(es) to consider for cluster analysis. Default: 'all'
    voxel_thresh : int, optional
        Threshold to apply to `stat_img`. If a negative number is provided a
        percentile threshold is used instead, where the percentile is
        determined by the equation `100 - voxel_thresh`. Default: 1.96
    cluster_extent : int, optional
        Minimum number of contiguous voxels required to consider a cluster in
        `filename`. Default: 20
    prob_thresh : [0, 100] int, optional
        Probability (percentage) threshold to apply to `atlas`, if it is
        probabilistic. Default: 5

    Returns
    -------
    clust_frame : pandas.DataFrame
        Dataframe wih information on clusters, including peak coordinates,
        average cluster value, volume of cluster, and percent overlap with
        neuroanatomical regions defined in `atlas`
    peaks_frame : pandas.DataFrame
        Dataframe with information on peaks, including peak coordinates, peak
        value, volume of cluster, and neuroanatomical location defined in
        `atlas`
    """
    stat_img = check_niimg(stat_img)
    atlas = check_atlases(atlas)
    voxel_volume = stat_img.header.get('pixdim')[1:4].prod()

    # threshold + clusterize image
    clust_img = process_img(stat_img,
                            voxel_thresh=voxel_thresh,
                            cluster_extent=cluster_extent)

    # get peak and cluster information
    clust_info, peaks_info = [], []
    for n, coord in enumerate(get_peak_coords(clust_img)):
        clust_data = image.index_img(clust_img, n).get_data()
        clust_mask = clust_data != 0
        peak_index = tuple(coord_xyz_to_ijk(clust_img.affine, coord))

        # basic info on cluster / peak
        peak_value = clust_data[peak_index]
        clust_mean = clust_data[clust_mask].mean()
        clust_volume = np.sum(clust_mask) * voxel_volume

        # atlas info on peak
        peak_info = get_peak_info(coord,
                                  atlastype=atlas,
                                  prob_thresh=prob_thresh)
        peak_summary = [
            peak if type(peak) != list else
            '; '.join(['{}% {}'.format(*e) for e in peak])
            for (_, peak) in peak_info
        ]

        # atlas info on cluster
        cluster_info = get_cluster_info(clust_mask,
                                        clust_img.affine,
                                        atlastype=atlas,
                                        prob_thresh=prob_thresh)
        clust_summary = [
            '; '.join(['{:.02f}% {}'.format(*e) for e in cluster])
            for (_, cluster) in cluster_info
        ]

        # append to output list
        clust_info += [coord + [clust_mean, clust_volume] + clust_summary]
        peaks_info += [coord + [peak_value, clust_volume] + peak_summary]

    clust_frame = pd.DataFrame(clust_info,
                               index=pd.Series(range(len(clust_info)),
                                               name='cluster_id'),
                               columns=['peak_x', 'peak_y', 'peak_z',
                                        'cluster_mean', 'volume'] + atlas)
    peaks_frame = pd.DataFrame(peaks_info,
                               index=pd.Series(range(len(peaks_info)),
                                               name='peak_id'),
                               columns=['peak_x', 'peak_y', 'peak_z',
                                        'peak_value', 'volume'] + atlas)

    return clust_frame, peaks_frame


def create_output(filename, atlas='all', voxel_thresh=1.96, cluster_extent=20,
                  prob_thresh=5, outdir=None):
    """
    Performs full cluster / peak analysis on `filename`

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
    atlas : str or list, optional
        Name of atlas(es) to consider for cluster analysis. Default: 'all'
    voxel_thresh : int, optional
        Threshold to apply to `stat_img`. If a negative number is provided a
        percentile threshold is used instead, where the percentile is
        determined by the equation `100 - voxel_thresh`. Default: 1.96
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

    # confirm input data is niimg_like to raise error as early as possible
    stat_img = check_niimg(filename)

    # get info for saving outputs
    if isinstance(filename, str):
        filename = op.abspath(filename)
        out_fname = op.basename(filename).split('.')[0]
        if outdir is None:
            outdir = op.dirname(filename)
    else:
        out_fname = 'mniatlasreader'
        if outdir is None:
            outdir = os.getcwd()

    # create output directory
    os.makedirs(outdir, exist_ok=True)

    # get cluster + peak information from image
    clust_frame, peaks_frame = get_statmap_info(stat_img, atlas=atlas,
                                                voxel_thresh=voxel_thresh,
                                                cluster_extent=cluster_extent,
                                                prob_thresh=prob_thresh)

    # write output .csv files
    clust_frame.to_csv(op.join(outdir, '{}_clusters.csv'.format(out_fname)))
    peaks_frame.to_csv(op.join(outdir, '{}_peaks.csv'.format(out_fname)))

    # generate stat map for plotting by collapsing all clusters into one image
    clust_img = process_img(stat_img,
                            voxel_thresh=voxel_thresh,
                            cluster_extent=cluster_extent)
    thresh_img = image.math_img('np.sum(img, axis=-1)', img=clust_img)

    # plot glass brain
    color_max = np.abs(thresh_img.get_data()).max()
    glass_fname = op.join(outdir, '{}.png'.format(out_fname))
    with warnings.catch_warnings():  # get rid of pesky warnings
        warnings.filterwarnings('ignore', category=FutureWarning)
        plotting.plot_glass_brain(thresh_img, vmax=color_max,
                                  threshold='auto', display_mode='lyrz',
                                  plot_abs=False, colorbar=True,
                                  black_bg=True,
                                  cmap=plotting.cm.cold_hot,
                                  output_file=glass_fname)

    # get template image for plotting cluster maps
    bgimg = nb.load(
        resource_filename(
            'mni_atlas_reader',
            'data/templates/MNI152_T1_1mm_brain.nii.gz'
        )
    )
    # plot clusters
    coords = clust_frame[['peak_x', 'peak_y', 'peak_z']].get_values()
    for idx, coord in enumerate(coords):
        clust_fname = '{}_cluster{:02d}.png'.format(out_fname, idx + 1)
        try:
            plotting.plot_stat_map(thresh_img, vmax=color_max,
                                   colorbar=True, title=clust_fname[:-4],
                                   threshold=voxel_thresh, draw_cross=True,
                                   black_bg=True, symmetric_cbar=True,
                                   output_file=op.join(outdir, clust_fname),
                                   bg_img=bgimg, cut_coords=coord,
                                   display_mode='ortho')
        except ValueError:
            plotting.plot_stat_map(thresh_img, vmax=color_max,
                                   colorbar=True, title=clust_fname[:-4],
                                   threshold=voxel_thresh, draw_cross=True,
                                   black_bg=True, symmetric_cbar=True,
                                   output_file=op.join(outdir, clust_fname))
