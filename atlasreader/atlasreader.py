"""
Primary functions of atlasreader
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
from scipy.ndimage import label, center_of_mass
from scipy.spatial.distance import cdist
from skimage.feature import peak_local_max
from sklearn.utils import Bunch


_ATLASES = [
    'aal',
    'aicha',
    'desikan_killiany',
    'destrieux',
    'harvard_oxford',
    'juelich',
    'marsatlas',
    'neuromorphometrics',
    'talairach_ba',
    'talairach_gyrus',
]

_DEFAULT = [
    'aal',
    'desikan_killiany',
    'harvard_oxford',
]


def get_atlas(atlastype, cache=True):
    """
    Gets `atlastype` image and corresponding label file from package resources

    Parameters
    ----------
    atlastype : str
        Name of atlas to query
    cache : bool, optional
        Whether to pre-load atlas image data. Default: True

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
    data_dir = resource_filename('atlasreader', 'data/atlases')
    atlas_path = op.join(data_dir, 'atlas_{0}.nii.gz'.format(atlastype))
    label_path = op.join(data_dir, 'labels_{0}.csv'.format(atlastype))

    if not all(op.exists(p) for p in [atlas_path, label_path]):
        raise ValueError('{} is not a valid atlas. Please check inputs and '
                         'try again.'.format(atlastype))

    atlas = Bunch(atlas=atlastype,
                  image=nb.load(atlas_path),
                  labels=pd.read_csv(label_path))
    if cache:
        atlas.image.get_data()

    return atlas


def check_atlases(atlases):
    """
    Checks atlases

    Parameters
    ----------
    atlases : str or list
        Name of atlas(es) to use

    Returns
    -------
    atlases : list
        Names of atlas(es) to use
    """
    if isinstance(atlases, str):
        atlases = [atlases]
    elif isinstance(atlases, dict):
        if all(hasattr(atlases, i) for i in ['image', 'atlas', 'labels']):
            return atlases
    if 'all' in atlases:
        draw = _ATLASES
    elif 'default' in atlases:
        draw = _DEFAULT
    else:
        draw = atlases

    return [get_atlas(a) if isinstance(a, str) else a for a in draw]


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
    labels = check_atlases(atlastype).labels
    try:
        return labels.query('index == {}'.format(label_id)).name.iloc[0]
    except IndexError:
        return 'no_label'


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
    xyz : (N, 3) numpy.ndarray
        Provided `coords` in `affine` space
    """
    coords = _check_coord_inputs(coords)
    mni_coords = np.dot(affine, coords)[:3].T
    return mni_coords


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
    ijk : (N, 3) numpy.ndarray
        Provided `coords` in cartesian space
    """
    coords = _check_coord_inputs(coords)
    vox_coords = np.linalg.solve(affine, coords)[:3].T.astype(int)
    return vox_coords


def get_peak_coords(clust_img):
    """
    Gets MNI coordinates of peak voxels within each cluster of `clust_img`

    Parameters
    ----------
    clust_img : 4D-niimg_like
        4D image of brain regions, where each volume is a separated cluster

    Returns
    ------
    coords : (N, 3) numpy.ndarray
        Coordinates of peak voxels in `clust_img`
    """

    # check cluster image and make it 4D, if not already
    clust_img = check_niimg(clust_img, atleast_4d=True)

    # create empty arrays to hold cluster size + peak coordinates
    clust_size = np.zeros(clust_img.shape[-1])
    maxcoords = np.zeros((clust_img.shape[-1], 3))

    # iterate through clusters and get info
    for n, cluster in enumerate(image.iter_img(clust_img)):
        cluster = np.abs(cluster.get_data())
        clust_size[n] = np.sum(cluster != 0)
        maxcoords[n] = center_of_mass(cluster == cluster.max())

    # sort peak coordinates by cluster size
    maxcoords = np.floor(maxcoords)[np.argsort(clust_size)[::-1]]

    # convert coordinates to MNI space
    coords = coord_ijk_to_xyz(clust_img.affine, maxcoords)

    return coords


def get_subpeak_coords(clust_img, min_distance=20):
    """
    Finds subpeaks in `clust_img` that are at least `min_distance` apart

    Parameters
    ----------
    clust_img : niimg_like
        Cluster image that `peaks` were derived from
    min_distance : float, optional
        Minimum distance required between peaks in `affine` (i.e., mm) space.
        Default: 20

    Returns
    -------
    peaks : (N, 3) numpy.ndarray
        Coordiantes of (sub)peak voxels in `clust_img`
    """
    data = check_niimg(clust_img).get_data()

    # find local maxima, excluding peaks that are on the border of the cluster
    local_max = peak_local_max(data, exclude_border=1, indices=False)

    # make new clusters to check for "flat" peaks + find CoM of those clusters
    labels, nl = label(local_max)
    ijk = center_of_mass(data, labels=labels, index=range(1, nl + 1))
    ijk = np.asarray(ijk, dtype=int)

    if len(ijk) > 1:
        # sort coordinates based on peak value
        ijk = ijk[data[tuple(map(tuple, ijk.T))].argsort()[::-1]]

        # convert to MNI space and get distance (in millimeters, not voxels)
        xyz = coord_ijk_to_xyz(clust_img.affine, ijk)
        distances = cdist(xyz, xyz)

        # remove "weaker" peak if it is too close to "stronger" peak
        keep = np.ones(len(xyz), dtype=bool)
        for r_idx, dist in enumerate(distances):
            if keep[r_idx] == 1:
                ind, = np.where(dist < min_distance)
                keep[np.setdiff1d(ind, r_idx)] = 0

        ijk = ijk[keep]

    coords = coord_ijk_to_xyz(clust_img.affine, ijk)

    return coords


def check_atlas_bounding_box(voxIDs, box_shape):
    """
    Returns the provided voxel ID if the voxel is inside the bounding box of
    the atlas image, otherwise the voxel ID will be replaced with the origin.

    Parameters
    ----------
    voxIDs : (N, 3) numpy.ndarray
        `coords` in cartesian space
    box_shape : (3,) list of int
        size of the atlas bounding box

    Returns
    ------
    ijk : (N, 3) numpy.ndarray
        `coords` in cartesian space that are inside the bounding box
    """

    # Detect voxels that are outside the atlas bounding box
    vox_outside_box = np.sum(
        (voxIDs < 0) + (voxIDs >= box_shape[:3]), axis=-1, dtype='bool')

    # Set those voxels to the origin (i.e. a voxel outside the brain)
    voxIDs[vox_outside_box] = np.zeros(3, dtype='int')

    return voxIDs


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

    # get atlas data
    atlastype = check_atlases(atlastype)
    data = atlastype.image.get_data()

    # get voxel index
    voxID = coord_xyz_to_ijk(atlastype.image.affine, coordinate).squeeze()
    voxID = check_atlas_bounding_box(voxID, data.shape)

    # get label information
    # probabilistic atlas is requested
    if atlastype.atlas.lower() in ['juelich', 'harvard_oxford']:
        probs = data[voxID[0], voxID[1], voxID[2]]
        probs[probs < prob_thresh] = 0
        idx, = np.where(probs)

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

    # get atlas data
    atlastype = check_atlases(atlastype)
    data = atlastype.image.get_data()

    # get coordinates of each voxel in cluster
    coords_vox = np.rollaxis(np.array(np.where(cluster)), 1)
    coords = coord_ijk_to_xyz(affine, coords_vox)

    # get voxel indexes
    voxIDs = coord_xyz_to_ijk(atlastype.image.affine, coords)
    voxIDs = check_atlas_bounding_box(voxIDs, data.shape)
    voxIDs = tuple(map(tuple, voxIDs.T))

    # get label information
    if atlastype.atlas.lower() in ['juelich', 'harvard_oxford']:
        labelIDs = np.argmax(data[voxIDs], axis=1)
        labelIDs[data[voxIDs].sum(-1) == 0] = -1
    else:
        labelIDs = data[voxIDs]

    unique_labels = np.unique(labelIDs)
    labels = np.array([get_label(atlastype, u) for u in unique_labels])
    N = float(len(labelIDs))
    percentage = np.array(
        [100 * np.sum(labelIDs == u) / N for u in unique_labels])

    sortID = np.argsort(percentage)[::-1]

    return [[percentage[s], labels[s]] for s in sortID if
            percentage[s] >= prob_thresh]


def process_img(stat_img, cluster_extent, voxel_thresh=1.96):
    """
    Parameters
    ----------
    stat_img : Niimg_like object
        Thresholded statistical map image
    cluster_extent : int
        Minimum number of voxels required to consider a cluster
    voxel_thresh : int, optional
        Threshold to apply to `stat_img`. If a negative number is provided a
        percentile threshold is used instead, where the percentile is
        determined by the equation `100 - voxel_thresh`. Default: 1.96

    Returns
    -------
    cluster_img : Nifti1Image
        4D image of brain regions, where each volume is a distinct cluster
    """
    # get input data image
    img_4d = check_niimg(stat_img, atleast_4d=True)
    if img_4d.shape[-1] == 1:
        stat_img = img_4d.slicer[..., 0]
    else:
        stat_img = image.index_img(img_4d, 0)

    # threshold image
    if voxel_thresh < 0:
        voxel_thresh = '{}%'.format(100 + voxel_thresh)
    else:
        # ensure that threshold is not greater than most extreme value in image
        if voxel_thresh > np.nan_to_num(np.abs(stat_img.get_data())).max():
            empty = np.zeros(stat_img.shape + (1,))
            return image.new_img_like(stat_img, empty)
    thresh_img = image.threshold_img(stat_img, threshold=voxel_thresh)

    # extract clusters
    min_region_size = cluster_extent * np.prod(thresh_img.header.get_zooms())
    clusters = []
    for sign in ['pos', 'neg']:
        # keep only data of given sign
        data = thresh_img.get_data().copy()
        data[(data < 0) if sign == 'pos' else (data > 0)] = 0

        # Do nothing if data array contains only zeros
        if np.any(data):
            try:
                if min_region_size != 0.0:
                    min_region_size -= 1e-8
                clusters += [connected_regions(
                    image.new_img_like(thresh_img, data),
                    min_region_size=min_region_size,
                    extract_type='connected_components')[0]]
            except TypeError:  # for no clusters
                pass

    # Return empty image if no clusters were found
    if len(clusters) == 0:
        return image.new_img_like(thresh_img, np.zeros(data.shape + (1,)))

    # Reorder clusters by their size
    clust_img = image.concat_imgs(clusters)
    cluster_size = (clust_img.get_data() != 0).sum(axis=(0, 1, 2))
    new_order = np.argsort(cluster_size)[::-1]
    clust_img_ordered = image.index_img(clust_img, new_order)

    return clust_img_ordered


def get_peak_data(clust_img, atlas='default', prob_thresh=5,
                  min_distance=None):
    """
    Parameters
    ----------
    clust_img : Niimg_like
        3D image of a single, valued cluster
    atlas : str or list, optional
        Name of atlas(es) to consider for cluster analysis. Default: 'default'
    prob_thresh : [0, 100] int, optional
        Probability (percentage) threshold to apply to `atlas`, if it is
        probabilistic. Default: 5
    min_distance : float, optional
        Specifies the minimum distance required between sub-peaks in a cluster.
        If None, sub-peaks will not be examined and only the primary cluster
        peak will be reported. Default: None

    Returns
    -------
    peak_summary : numpy.ndarray
        Info on peaks found in `clust_img`, including peak coordinates, peak
        values, volume of cluster that peaks belong to, and neuroanatomical
        locations defined in `atlas`
    """
    data = check_niimg(clust_img).get_data()

    # get voxel volume information
    voxel_volume = np.prod(clust_img.header.get_zooms())

    # find peaks -- or subpeaks, if `min_distance` is set
    if min_distance is None:
        peaks = get_peak_coords(clust_img)
    else:
        peak_img = image.math_img('np.abs(img)', img=clust_img)
        peaks = get_subpeak_coords(peak_img, min_distance=min_distance)
    peaks = coord_xyz_to_ijk(clust_img.affine, peaks)

    # get info on peak in cluster (all peaks belong to same cluster ID!)
    cluster_volume = np.repeat(np.sum(data != 0) * voxel_volume, len(peaks))
    peak_values = data[tuple(map(tuple, peaks.T))]
    coords = coord_ijk_to_xyz(clust_img.affine, peaks)
    peak_info = []
    for coord in coords:
        coord_info = []
        for atype in check_atlases(atlas):
            segment = read_atlas_peak(atype, coord, prob_thresh)
            coord_info.append([atype.atlas, segment])

        peak_info += [[peak if type(peak) != list else
                       '; '.join(['{}% {}'.format(*e) for e in peak])
                       for (_, peak) in coord_info]]

    return np.column_stack([coords, peak_values, cluster_volume, peak_info])


def get_cluster_data(clust_img, atlas='default', prob_thresh=5):
    """
    Parameters
    ----------
    clust_img : Niimg_like
        3D image of a single, valued cluster
    atlas : str or list, optional
        Name of atlas(es) to consider for cluster analysis. Default: 'default'
    prob_thresh : [0, 100] int, optional
        Probability (percentage) threshold to apply to `atlas`, if it is
        probabilistic. Default: 5

    Returns
    -------
    clust_frame : list
        Info on cluster in `clust_img`, including coordinates of peak,
        average cluster value, volume of cluster, and percent overlap with
        neuroanatomical regions defined in `atlas`
    """
    data = check_niimg(clust_img).get_data()
    voxel_volume = np.prod(clust_img.header.get_zooms())

    coord = get_peak_coords(clust_img).squeeze().tolist()
    clust_mask = data != 0
    clust_mean = data[clust_mask].mean()
    cluster_volume = np.sum(clust_mask) * voxel_volume

    # atlas info on cluster
    cluster_info = []
    for atype in check_atlases(atlas):
        segment = read_atlas_cluster(atype, clust_mask,
                                     clust_img.affine, prob_thresh)
        cluster_info.append([atype, segment])

    cluster_info = ['; '.join(['{:.02f}% {}'.format(*e) for e in cluster])
                    for (_, cluster) in cluster_info]

    return coord + [clust_mean, cluster_volume] + cluster_info


def get_statmap_info(stat_img, cluster_extent, atlas='default',
                     voxel_thresh=1.96, prob_thresh=5, min_distance=None):
    """
    Extract peaks and cluster information from `clust_img` for `atlas`

    Parameters
    ----------
    stat_img : Niimg_like
        4D image of brain regions, where each volume is a distinct cluster
    cluster_extent : int
        Minimum number of contiguous voxels required to consider a cluster in
        `stat_img`
    atlas : str or list, optional
        Name of atlas(es) to consider for cluster analysis. Default: 'default'
    voxel_thresh : int, optional
        Threshold to apply to `stat_img`. If a negative number is provided a
        percentile threshold is used instead, where the percentile is
        determined by the equation `100 - voxel_thresh`. Default: 1.96
    prob_thresh : [0, 100] int, optional
        Probability (percentage) threshold to apply to `atlas`, if it is
        probabilistic. Default: 5
    min_distance : float, optional
        Specifies the minimum distance required between sub-peaks in a cluster.
        If None, sub-peaks will not be examined and only the primary cluster
        peak will be reported. Default: None

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
    atlas = check_atlases(atlas)  # loading in the atlases takes the longest

    # threshold + clusterize image
    clust_img = process_img(stat_img,
                            voxel_thresh=voxel_thresh,
                            cluster_extent=cluster_extent)

    clust_info, peaks_info = [], []
    if clust_img.get_data().any():
        for n, cluster in enumerate(image.iter_img(clust_img)):
            peak_data = get_peak_data(cluster, atlas=atlas,
                                      prob_thresh=prob_thresh,
                                      min_distance=min_distance)
            clust_data = get_cluster_data(cluster, atlas=atlas,
                                          prob_thresh=prob_thresh)

            cluster_id = np.repeat(n + 1, len(peak_data))
            peaks_info += [np.column_stack([cluster_id, peak_data])]
            clust_info += [[n + 1] + clust_data]
        clust_info = np.row_stack(clust_info)
        peaks_info = np.row_stack(peaks_info)

    # construct dataframes and reset floats
    atlasnames = [a.atlas for a in atlas]
    clust_frame = pd.DataFrame(clust_info,
                               columns=['cluster_id',
                                        'peak_x', 'peak_y', 'peak_z',
                                        'cluster_mean', 'volume_mm']
                               + atlasnames)
    peaks_frame = pd.DataFrame(peaks_info,
                               columns=['cluster_id',
                                        'peak_x', 'peak_y', 'peak_z',
                                        'peak_value', 'volume_mm']
                               + atlasnames)
    for col in range(6):
        clust_frame.iloc[:, col] = clust_frame.iloc[:, col].astype(float)
        peaks_frame.iloc[:, col] = peaks_frame.iloc[:, col].astype(float)

    return clust_frame, peaks_frame


def create_output(filename, cluster_extent, atlas='default', voxel_thresh=1.96,
                  prob_thresh=5, min_distance=None, outdir=None,
                  glass_plot_kws=None, stat_plot_kws=None):
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
    filename : Niimg_like
        A 3D statistical image.
    cluster_extent : int
        Minimum number of contiguous voxels required to consider a cluster in
        `filename`
    atlas : str or list, optional
        Name of atlas(es) to consider for cluster analysis. Default: 'default'
    voxel_thresh : int, optional
        Threshold to apply to `stat_img`. If a negative number is provided a
        percentile threshold is used instead, where the percentile is
        determined by the equation `100 - voxel_thresh`. Default: 1.96
    prob_thresh : int, optional
        Probability (percentage) threshold to apply to `atlas`, if it is
        probabilistic. Default: 5
    min_distance : float, optional
        Specifies the minimum distance (in mm) required between sub-peaks in a
        cluster. If None, sub-peaks will not be examined and only the primary
        cluster peak will be reported. Default: None
    outdir : str or None, optional
        Path to desired output directory. If None, generated files will be
        saved to the same folder as `filename`. Default: None
    glass_plot_kws : dict or None, optional
        Additional keyword arguments to pass to
        `nilearn.plotting.plot_glass_brain`.
    stat_plot_kws : dict or None, optional
        Additional keyword arguments to pass to
        `nilearn.plotting.plot_stat_map`.
    """

    # confirm input data is niimg_like to raise error as early as possible
    stat_img = check_niimg(filename)

    # get info for saving outputs
    if isinstance(filename, str):
        filename = op.abspath(filename)
        if filename.endswith('.nii.gz'):
            out_fname = op.basename(filename)[:-7]
        elif filename.endswith('.nii'):
            out_fname = op.basename(filename)[:-4]
        elif filename.endswith('.img'):
            out_fname = op.basename(filename)[:-4]
        if outdir is None:
            outdir = op.dirname(filename)
    else:
        out_fname = 'atlasreader'
        if outdir is None:
            outdir = os.getcwd()

    # create output directory
    os.makedirs(outdir, exist_ok=True)

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

        glass_plot_params = {
            'stat_map_img': thresh_img,
            'output_file': glass_fname,
            'display_mode': 'lyrz',
            'colorbar': True,
            'black_bg': True,
            'cmap': plotting.cm.cold_hot,
            'vmax': color_max,
            'plot_abs': False,
            'symmetric_cbar': False
        }
        if glass_plot_kws is None:
            glass_plot_kws = {}
        glass_plot_params.update(glass_plot_kws)
        plotting.plot_glass_brain(**glass_plot_params)

    # Check if thresholded image contains only zeros
    if np.any(thresh_img.get_data()):

        # get cluster + peak information from image
        clust_frame, peaks_frame = get_statmap_info(
            stat_img, atlas=atlas, voxel_thresh=voxel_thresh,
            cluster_extent=cluster_extent, prob_thresh=prob_thresh,
            min_distance=min_distance)

        # write output .csv files
        clust_frame.to_csv(op.join(
            outdir, '{}_clusters.csv'.format(out_fname)),
            index=False, float_format='%5g')
        peaks_frame.to_csv(op.join(
            outdir, '{}_peaks.csv'.format(out_fname)),
            index=False, float_format='%5g')

        # get template image for plotting cluster maps
        bgimg = nb.load(
            resource_filename(
                'atlasreader',
                'data/templates/MNI152_T1_1mm_brain.nii.gz'
            )
        )
        # plot clusters
        coords = clust_frame[['peak_x', 'peak_y', 'peak_z']].get_values()
        for idx, coord in enumerate(coords):
            clust_fname = '{}_cluster{:02d}.png'.format(out_fname, idx + 1)
            stat_plot_params = {
                'stat_map_img': thresh_img,
                'bg_img': bgimg,
                'cut_coords': coord,
                'output_file': op.join(outdir, clust_fname),
                'colorbar': True,
                'title': clust_fname[:-4],
                'threshold': voxel_thresh,
                'black_bg': True,
                'symmetric_cbar': False,
                'vmax': color_max
            }
            if stat_plot_kws is None:
                stat_plot_kws = {}
            stat_plot_params.update(stat_plot_kws)
            plotting.plot_stat_map(**stat_plot_params)
