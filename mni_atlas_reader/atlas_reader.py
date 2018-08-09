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
    'Desikan-Killiany', 'Juelich', 'Neuromorphometrics', 'Destrieux',
    'Harvard_Oxford', 'AAL'
]
_ACCEPTED_ATLASES = _ATLASES + [a.lower() for a in _ATLASES]


def get_atlas(atlas):
    """
    Gets `atlas` image and label file from package resources

    Parameters
    ----------
    atlas : str
        Name of atlas for which to get information

    Returns
    -------
    info : sklearn.utils.Bunch
        image : Niimg_like
            ROI image loaded and ready for use
        labels : pandas.core.data.DataFrame
            Dataframe with columns ['index', 'name'] matching regions in
            `image` to anatomical names
    """
    # we accept uppercase atlases via argparse but our filenames are lowercase
    atlas = atlas.lower()

    # get the path to atlas + label files shipped with package
    # resource_filename ensures that we're getting the correct path
    data_dir = resource_filename('mni_atlas_reader', 'data/atlases')
    atlas_path = op.join(data_dir, 'atlas_{0}.nii.gz'.format(atlas))
    label_path = op.join(data_dir, 'labels_{0}.csv'.format(atlas))

    # return loaded filenames (we should only have to call this once!)
    return Bunch(image=nb.load(atlas_path), labels=pd.read_csv(label_path))


def get_vox_coord(affine, coord):
    """
    Computes voxel matrix index from MNI coordinate.

    Parameters
    ----------
    affine : (4, 4) array-like
        Affine matrix.
    coord : list of floats
        x, y, z MNI coordinates

    Returns
    ------
    coordinates : list
        x, y, z array index of voxel
    """
    inverse = np.linalg.inv(affine)
    vox_coord = np.dot(inverse, np.hstack((coord, 1)))[:3]
    return vox_coord.round().astype('int').tolist()


def get_label(atlas, label_id):
    """Reads out the name of a specific label

    Parameters
    ----------
    atlas : str
        Name of atlas to use
    label_id : int
        Numeric label of cluster

    Returns
    ------
    label : str
        Atlas label name
    """

    labels = get_atlas(atlas).labels
    try:
        return labels.query('index=={}'.format(label_id)).name.iloc[0]
    except IndexError:
        return 'No_label'


def get_clusters(data, min_extent=5):
    """
    Extracts out clusters from a stat map.

    Parameters
    ----------
    data : (X, Y, Z) array-like
        Stat map array
    min_extent : int
        Minimum extension of cluster

    Returns
    ------
    clusters : ndarray
        Array of numerically labelled clusters
    nclusters : int
        Number of clusters
    """
    clusters, nclusters = label(data)
    for idx in range(1, nclusters + 1):
        if np.sum(clusters == idx) < min_extent:
            clusters[clusters == idx] = 0
    nclusters = len(np.setdiff1d(np.unique(clusters), [0]))
    return clusters, nclusters


def get_peak_coords(img, affine, data):
    """
    Gets MNI coordinates of peak voxel within each cluster.

    Parameters
    ----------
    img : (X, Y, Z) array-like
        Array of numerically labelled clusters
    affine : (4, 4) array-like
        NIfTI affine matrix
    data : (X, Y, Z) array-like
        Stat map array

    Returns
    ------
    coordinates : list of lists
        x, y, z MNI coordinates of peak voxels
    """
    coords = []
    clusters = np.setdiff1d(np.unique(img.ravel()), [0])
    cs = []
    maxcoords = []
    for lab in clusters:
        cs.append(np.sum(img == lab))
        maxval = np.max(data[img == lab])
        maxidx = np.nonzero(np.multiply(data, img == lab) == maxval)
        maxcoords.append([m[0] for m in maxidx])

    maxcoords = np.asarray(maxcoords)
    maxcoords = maxcoords[np.argsort(cs)[::-1], :]
    for i, lab in enumerate(clusters[np.argsort(cs)[::-1]]):
        coords.append(np.dot(affine,
                             np.hstack((maxcoords[i], 1)))[:3].tolist())
    return coords


def get_cluster_coords(cluster, affine):
    """
    Get matrix indices of voxels in cluster.

    Parameters
    ----------
    cluster : (X, Y, Z) array-like
        Boolean mask of cluster
    affine : (4, 4) array-like
        NIfTI affine matrix

    Returns
    ------
    coordinates : list of numpy.ndarray
        List of coordinates for voxels in cluster
    """
    coords_vox = np.rollaxis(np.array(np.where(cluster)), 1)
    coords = [np.dot(affine, np.hstack((c, 1)))[:3] for c in coords_vox]
    return coords


def read_atlas_peak(atlastype, coordinate, probThresh=5):
    """
    Reads specific atlas and returns segment/probability information.

    It is possible to threshold a given probability atlas [in percentage].

    Parameters
    ----------
    atlastype : str
        Name of atlas to use
    coordinate : list of float
        x, y, z MNI coordinates of voxel
    probThresh : int
        Probability threshold for when using a probabilistic atlas

    Returns
    ------
    label : list or str
        Atlas label of peak coordinate. If using a probabilistic atlas,
        then returns a list of the probability information organized as
        [probability, label]. Otherwise, the atlas label is returned.
    """

    atlas = get_atlas(atlastype).image

    # Get atlas data and affine matrix
    data = atlas.get_data()
    affine = atlas.affine

    # Get voxel index
    voxID = get_vox_coord(affine, coordinate)

    # Get Label information
    if atlastype in ['Juelich', 'Harvard_Oxford']:
        probs = data[voxID[0], voxID[1], voxID[2]]
        probs[probs < probThresh] = 0
        idx = np.where(probs)[0]

        # sort list by probability
        idx = idx[np.argsort(probs[idx])][::-1]

        # get probability and label names
        probLabel = []
        for i in idx:
            label = get_label(atlastype, i)
            probLabel.append([probs[i], label])

        # If no labels found
        if probLabel == []:
            probLabel = [[0, 'No_label']]

        return probLabel

    else:
        labelID = int(data[voxID[0], voxID[1], voxID[2]])
        label = get_label(atlastype, labelID)
        return label


def read_atlas_cluster(atlastype, cluster, affine, probThresh):
    """
    Reads specific atlas and returns segment/probability information.

    It is possible to threshold a given probability atlas (in percentage).

    Parameters
    ----------
    atlastype : str
        Name of atlas to use.
    cluster : (X, Y, Z) array-like
        Boolean mask of cluster
    affine : (4, 4) array-like
        NIfTI affine matrix
    probThresh : int
        Probability threshold for when using a probabilistic atlas

    Returns
    ------
    segments : list of lists
        List of probability information for the cluster, which is
        organized as [probability, region name]
    """
    atlas = get_atlas(atlastype).image

    # Get atlas data and affine matrix
    atlas_data = atlas.get_data()
    atlas_affine = atlas.affine

    # Get coordinates of each voxel in cluster
    coords = get_cluster_coords(cluster, affine)

    # Get voxel indexes
    voxIDs = [get_vox_coord(atlas_affine, c) for c in coords]

    # Get Label information
    if atlastype in ['Juelich', 'Harvard_Oxford']:
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
            percentage[s] >= probThresh]


def get_peak_info(coord, atlastype='all', probThresh=5):
    """
    Gets region and probability information for a peak voxel based on the
    provided atlas.

    Parameters
    ----------
    coord : list of floats
        x, y, z MNI coordinates of voxel
    atlastype : str
        Name of atlas to use
    probThresh : int
        Probability threshold for when using a probabilistic atlas

    Returns
    ------
    peakinfo :
        List of lists containing atlas name and probability information
        for the atlas labels associated with the voxel
    """
    peakinfo = []
    if isinstance(atlastype, str):
        atlastype = [atlastype]

    if 'all' not in atlastype:
        for atypes in atlastype:
            segment = read_atlas_peak(atypes, coord, probThresh)
            peakinfo.append([atypes, segment])
    else:
        for atypes in _ATLASES:
            segment = read_atlas_peak(atypes, coord, probThresh)
            peakinfo.append([atypes, segment])

    return peakinfo


def get_cluster_info(cluster, affine, atlastype='all', probThresh=5):
    """
    Gets region and probability information for a cluster based on the
    provided atlas.

    Parameters
    ----------
    cluster :
        Cluster number
    affine : (4, 4) array-like
        NIfTI affine matrix
    atlastype : str
        Atlas name
    probThresh : int
        Probability threshold for when using a probabilistic atlas

    Returns
    ------
    clusterinfo : list of lists
        List containing the cluster information, which includes the
        atlas name, and the probability information of each atlas label
        for cluster.
    """

    clusterinfo = []

    if isinstance(atlastype, str):
        atlastype = [atlastype]

    if 'all' not in atlastype:
        for atypes in atlastype:
            segment = read_atlas_cluster(atypes, cluster, affine, probThresh)
            clusterinfo.append([atypes, segment])
    else:
        for atypes in _ATLASES:
            segment = read_atlas_cluster(atypes, cluster, affine, probThresh)
            clusterinfo.append([atypes, segment])

    return clusterinfo


def create_output(filename, atlas, voxelThresh=2, clusterExtend=5,
                  probabilityThreshold=5, outDir=None):
    """
    Generates output table containing each clusters' number of voxels,
    average activation across voxels, peak voxel coordinates, and
    neuroanatomical location of the peak voxel based on the specified atlas.

    In addition, separate stat maps are created for each cluster. In each
    image, cross hairs are located on the cluster's peak voxel.

    Parameters
    ----------
    filename : str
        The full or relative path to the statistical map to use
    atlas : str
        Atlas name to use
    voxelThresh : int
        Value threshold for voxels to be considered in cluster extraction
    clusterExtend : int
        Required number of contiguous voxels for a cluster to be retained for
        analysis
    probabilityThreshold : int
        Probability threshold for when using a probabilistic atlas
    outDir : str or None
        Path to the output directory. If None is specified, output will be
        stored in the same folder as the input files.

    Returns
    ------
    None
    """

    fname = os.path.abspath(filename)

    # set up output directory
    if outDir is not None:
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        savedir = outDir

    else:
        savedir = os.path.dirname(fname)

    # set up output file name which is the same name (no extension) of input
    # file
    out_fname = os.path.basename(fname).split('.')[0]

    # Get data from NIfTI file
    img = nb.load(filename)
    imgdata = img.get_data()
    if len(imgdata.shape) != 3:
        imgdata = imgdata[:, :, :, 0]

    # Get top x-% of voxels if voxelThresh is negative
    if voxelThresh < 0:
        voxelThresh = np.percentile(
            np.abs(imgdata[imgdata != 0]), (100 + voxelThresh))

    # Get clusters from data
    clusters, nclusters = get_clusters(
        abs(imgdata) > voxelThresh, min_extent=clusterExtend)

    # Clean img data
    imgdata[clusters == 0] = 0
    new_image = nb.Nifti1Image(imgdata, img.affine, img.header)

    # Plot Glass Brain
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

    # Get coordinates of peaks
    coords = get_peak_coords(clusters, img.affine, np.abs(imgdata))

    # Get Peak and Cluster information
    peak_summary = []
    peak_value = []
    cluster_summary = []
    cluster_mean = []
    volume_summary = []

    for n, coord in enumerate(coords):
        peakinfo = get_peak_info(
            coord, atlastype=atlas, probThresh=probabilityThreshold)
        peak_summary.append([p[1] if type(p[1]) != list else '; '.join(
            ['% '.join(e) for e in np.array(p[1])]) for p in peakinfo])
        voxID = get_vox_coord(img.affine, coord)
        peak_value.append(imgdata[voxID[0], voxID[1], voxID[2]])

        idx = get_vox_coord(img.affine, coord)
        clusterID = clusters[idx[0], idx[1], idx[2]]
        clusterinfo = get_cluster_info(
            clusters == clusterID,
            img.affine,
            atlastype=atlas,
            probThresh=probabilityThreshold)
        cluster_summary.append(['; '.join(
            ['% '.join([str(round(e[0], 2)), e[1]]) for e in coord[1]]) for
             coord in clusterinfo])
        cluster_mean.append(imgdata[clusters == clusterID].mean())

        voxel_volume = int(img.header['pixdim'][1:4].prod())
        volume_summary.append(np.sum(clusters == clusterID) * voxel_volume)

    # Write output .csv file
    header = [p[0] for p in peakinfo]
    with open(os.path.join(savedir, '{}.csv'.format(out_fname)), 'w') as f:
        f.writelines(','.join(
            ['ClusterID', 'Peak_Location', 'Cluster_Mean', 'Volume'] + header)
            + '\n')

        for i, c in enumerate(cluster_summary):
            f.writelines(
                ','.join(['Cluster%.02d' % (i + 1), '_'.join(
                    [str(xyz) for xyz in coords[i]]), str(cluster_mean[i]),
                     str(volume_summary[i])] + c) + '\n')

        f.writelines('\n')

        f.writelines(
            ','.join(['PeakID', 'Peak_Location', 'Peak_Value', 'Volume'] +
                     header) + '\n')

        for i, p in enumerate(peak_summary):
            f.writelines(
                ','.join(['Peak%.02d' % (i + 1), '_'.join(
                    [str(xyz) for xyz in coords[i]]), str(peak_value[i]),
                     str(volume_summary[i])] + p) + '\n')

    # get template image for plotting cluster maps
    bgimg = nb.load(
        resource_filename(
            'mni_atlas_reader',
            'data/templates/MNI152_T1_1mm_brain.nii.gz'
        )
    )
    # Plot Clusters
    for idx, coord in enumerate(coords):
        cluster_name = '{}_cluster0{}'.format(out_fname, idx + 1)
        out_cluster_file = os.path.join(savedir, '{}.png'.format(cluster_name))

        try:
            plot_stat_map(new_image, bg_img=bgimg, cut_coords=coord,
                          display_mode='ortho', colorbar=True,
                          title=cluster_name, threshold=voxelThresh,
                          draw_cross=True, black_bg=True, symmetric_cbar=True,
                          vmax=color_max, output_file=out_cluster_file)
        except ValueError:
            plot_stat_map(new_image, vmax=color_max,
                          colorbar=True, title=cluster_name,
                          threshold=voxelThresh, draw_cross=True,
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
    """ Reads command line arguments and returns """
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
                        dest='voxelThresh', metavar='threshold',
                        help='Value threshold that voxels in provided file '
                             'must surpass in order to be considered in '
                             'cluster extraction.')
    parser.add_argument('-c', '--cluster', type=float, default=5,
                        dest='clusterExtend', metavar='extent',
                        help='Required number of contiguous voxels for a '
                             'cluster to be retained for analysis.')
    parser.add_argument('-p', '--probability', type=check_limit, default=5,
                        dest='probabilityThreshold', metavar='threshold',
                        help='Threshold to consider when using a '
                             'probabilistic atlas for extracting anatomical '
                             'cluster locations. Value will apply to all '
                             'request probabilistic atlases, and should range '
                             'between 0 and 100.')
    parser.add_argument('-o', '--outdir', type=str, dest='outDir',
                        metavar='outdir',
                        help='Output directory for created files. If it is '
                        'not specified, then output files are created in the '
                        'same directory as the statistical map that is '
                        'provided.')

    return parser.parse_args()


def main():
    """
    The primary entrypoint for calling atlas reader via the command line

    All parameters are read via argparse, so this should only be called from
    the command line!
    """

    opts = _get_parser()
    create_output(opts.filename, opts.atlas,
                  voxelThresh=opts.voxelThresh,
                  clusterExtend=opts.clusterExtend,
                  probabilityThreshold=opts.probabilityThreshold,
                  outDir=opts.outDir)


if __name__ == '__main__':
    main()
