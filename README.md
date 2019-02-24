[![codecov](https://codecov.io/gh/miykael/atlasreader/branch/master/graph/badge.svg)](https://codecov.io/gh/miykael/atlasreader)
[![Build Status](https://travis-ci.org/miykael/atlasreader.svg?branch=master)](https://travis-ci.org/miykael/atlasreader)
[![GitHub issues](https://img.shields.io/github/issues/miykael/atlasreader.svg)](https://github.com/miykael/atlasreader/issues/)
[![GitHub pull-requests](https://img.shields.io/github/issues-pr/miykael/atlasreader.svg)](https://github.com/miykael/atlasreader/pulls/)
[![GitHub contributors](https://img.shields.io/github/contributors/miykael/atlasreader.svg)](https://GitHub.com/miykael/atlasreader/graphs/contributors/)
[![GitHub Commits](https://github-basic-badges.herokuapp.com/commits/miykael/atlasreader.svg)](https://github.com/miykael/atlasreader/commits/master)
[![GitHub size](https://github-size-badge.herokuapp.com/miykael/atlasreader.svg)](https://github.com/miykael/atlasreader/archive/master.zip)
[![GitHub HitCount](http://hits.dwyl.io/miykael/atlasreader.svg)](http://hits.dwyl.io/miykael/atlasreader)

# AtlasReader

This package provides a Python interface for generating coordinate tables and
region labels from statistical MRI images. It is intended for neuroscience
researchers and neuroimaging enthusiasts who are looking for a quick and easy way
to localize and extract relevant peak and cluster information and create
informative and nice looking overview figures.

Please check out our interactive notebook on mybinder.org to see `atlasreader`
in action:  
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/miykael/atlasreader/master?filepath=notebooks%2Fatlasreader.ipynb)

If you are using `atlasreader` in your publication, please cite the following paper:  
[![DOI](http://joss.theoj.org/papers/10.21105/joss.01257/status.svg)](https://doi.org/10.21105/joss.01257)

Notter M. P., Gale D., Herholz P., Markello R. D., Notter-Bielser M.-L., & Whitaker K. (2019). AtlasReader: A Python package to generate coordinate tables, region labels, and informative figures from statistical MRI images. *Journal of Open Source Software, 4(34), 1257*, [https://doi.org/10.21105/joss.01257](https://doi.org/10.21105/joss.01257).


## Installation

This package requires Python >= 3.6. Provided you have `pip` at your disposal,
installing `atlasreader` is as simple as this:

```bash
pip install atlasreader
```

If you want to build `atlasreader` directly from source code, use the 
following code:

```bash
git clone https://github.com/miykael/atlasreader.git
cd atlasreader
python setup.py install
```


## Usage

AtlasReader can either be run through the command line interface or directly
within Python. The commands to do so are rather straight forward. Let's say you
want to apply AtlasReader to a statistical image called
`file_name = 'stat_img.nii'`, and only want to keep clusters if they have more
than 5 voxels:

#### Python
```python
from atlasreader import create_output
create_output(file_name, cluster_extent=5)
```

#### Command Line
```bash
atlasreader file_name 5
```


### Outputs

After executing AtlasReader on a given image, four kinds of outputs are generated:

1. An **overview figure** that shows the results within the whole brain at once  
   ![Overview Figure](paper/fig_overview_figure.png)

2. **For each cluster**, an **informative figure** showing the sagittal, coronal
   and transversal plane centered on the main peak of the cluster  
   ![Cluster Figure](paper/fig_cluster_figure.png)

3. A **csv file** containing relevant information about the **peak** of each
   cluster. This table contains the cluster association and location of each
   peak, its signal value at this location, the cluster extent (in mm, not in
   number of voxels), as well as the membership of each peak, given a
   particular atlas.  
   ![Table Peak](paper/table_peak.png)

4. A **csv** file containing relevant information about each **cluster**. Table
   showing relevant information for the cluster extent of each ROI. This table
   contains the cluster association and location of each peak, the mean value
   within the cluster, the cluster extent (in mm, not in number of voxels), as
   well as the membership of each cluster, given a particular atlas.  
   ![Table Cluster](paper/table_cluster.png)


### Additional parameters

`atlasreader.create_output` has many additional parameters that allow you to change the way
the clusters are generated and what kind of outputs are generated:

- **filename**: Niimg_like  
    A 3D statistical image.
- **cluster_extent**: int  
    Minimum number of contiguous voxels required to consider a cluster in `filename`
- **atlas**: str or list, optional  
    Name of atlas(es) to consider for cluster analysis. ***Default***: `'default'`
- **voxel_thresh**: float, optional
    Threshold to apply to `stat_img`.  Use `direction` to specify the 
    directionality of the threshold. If a negative number is provided a
    percentile threshold is used instead, where the percentile is determined
    by the equation `100 - voxel_thresh`. ***Default***: `1.96`
- **direction**: str, optional
    Specifies the direction in which `voxel_thresh` should be applied. Possible
    values are `'both'`, `'pos'` or `'neg'`. ***Default***: `'both'`
- **prob_thresh**: int, optional  
    Probability (percentage) threshold to apply to `atlas`, if it is
    probabilistic. ***Default***: `5`
- **min_distance**: float, optional  
    Specifies the minimum distance (in mm) required between sub-peaks in a
    cluster. If None, sub-peaks will not be examined and only the primary
    cluster peak will be reported. ***Default***: `None`
- **outdir**: str or None, optional  
    Path to desired output directory. If None, generated files will be
    saved to the same folder as `filename`. ***Default***: `None`
- **glass_plot_kws**: dict or None, optional  
    Additional keyword arguments to pass to `nilearn.plotting.plot_glass_brain`.
    ***Default***: `None`
- **stat_plot_kws**: dict or None, optional  
    Additional keyword arguments to pass to `nilearn.plotting.plot_stat_map`.
    ***Default***: `None`

For a more detailed explanation about the toolbox and the effect of the
parameters above, see the [example notebook](https://github.com/miykael/atlasreader/blob/master/notebooks/atlasreader.ipynb).
You can checkout the notebook either interactively via [mybinder.org](https://mybinder.org/v2/gh/miykael/atlasreader/master?filepath=notebooks%2Fatlasreader.ipynb) or explore a read-only version on [nbviewer.jupyter.org](https://nbviewer.jupyter.org/github/miykael/atlasreader/blob/master/notebooks/atlasreader.ipynb).


## How to get involved

We're thrilled to welcome new contributors!

If you're interested in getting involved, you should start by reading our
[contributing guidelines](CONTRIBUTING.md).

Once you're done with that, you can take a look at our list of active
[issues](https://github.com/miykael/atlasreader/issues) and let us know if
there's something you'd like to begin working on.

If you've found a bug, are experiencing a problem, or have a question, create a
new [issue](https://github.com/miykael/atlasreader/issues) with some information
about it!


## Licence

AtlasReader is licensed under the BSD-3 license; however, the atlases it uses 
are separately licensed under more restrictive frameworks.
By using AtlasReader, you agree to abide by the license terms of the
individual atlases. Information on these terms can be found online at:
https://github.com/miykael/atlasreader/tree/master/atlasreader/data
