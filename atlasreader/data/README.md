# Template Information

## MNI152 T1 1mm Template

The MNI152 T1 template is taken from FSL 5.0. It's kindly supplied by Andrew
Janke. This template is derived from 152 structural images, averaged together
after high-dimensional nonlinear registration into the common MNI152 co-ordinate
system. It corresponds to the "[152 nonlinear 6th generation](http://www.bic.mni.mcgill.ca/ServicesAtlases/HomePage)" atlas.

### Creation of template

This atlas is a direct copy of the file `MNI152_T1_1mm_brain.nii.gz` packaged
with FSL 5.0.

### License

The FSL templates and atlases are under the License that can be found
[here](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Licence).

## ICBM 2009c Nonlinear Asymmetric Template

A number of unbiased non-linear averages of the MNI152 database have been
generated that combines the attractions of both high-spatial resolution and
signal-to-noise while not being subject to the vagaries of any single brain
(Fonov et al., 2011). The procedure involved multiple iterations of a process
where, at each iteration, individual native MRIs were non-linearly fitted to
the average template from the previous iteration, beginning with the MNI152
linear template. The templates present an unbiased standard magnetic resonance
imaging template brain volume for normal population. These volumes were created
using data from ICBM project. For more information see: http://www.bic.mni.mcgill.ca/ServicesAtlases/ICBM152NLin2009.

### Creation of template

The template was download on 24. Sept. 2018 from
http://www.bic.mni.mcgill.ca/~vfonov/icbm/2009/mni_icbm152_nlin_asym_09c_nifti.zip.
We used the following command to create a T1 brain only template:

```
fslmaths mni_icbm152_t1_tal_nlin_asym_09c.nii \
    -mul mni_icbm152_t1_tal_nlin_asym_09c_mask.nii \
    mni_icbm152_t1_tal_nlin_asym_09c_brain.nii.gz \
    -odt short
```

### License

Copyright (C) 1993-2004 Louis Collins, McConnell Brain Imaging Centre, Montreal
Neurological Institute, McGill University. Permission to use, copy, modify, and
distribute this software and its documentation for any purpose and without fee
is hereby granted, provided that the above copyright notice appear in all
copies. The authors and McGill University make no representations about the
suitability of this software for any purpose.  It is provided "as is" without
express or implied warranty.  The authors are not responsible for any data loss,
equipment damage, property loss, or injury to subjects or patients resulting
from the use or misuse of this software package.

# Atlas Information

`atlasreader` contains many different atlases that each are under their own
license and related to specific publications. This README file acknoledges those
licenses, describes the origin of the atlases and how they were acquired and
explains how the atlases were adapted to fit into the `atlasreader` framework.
Many of those atlases were downloaded with `nilearn.datasets`'s `fetch_`
functions. For more about this, go to [`nilearn`'s official homepage](http://nilearn.github.io/modules/reference.html#module-nilearn.datasets).

## Anatomical Automatic Labeling 2 (AAL2)

Anatomical Automatic Labeling (AAL) is a package for the anatomical labeling of
functional brain mapping experiments. It is an in-house package made by
Neurofunctional Imaging Group (GIN, UMR5296, Bordeaux, France), which is
available to the scientific community as a copyright freeware under the terms of
the GNU General Public License. For more information on this dataset’s
structure, see http://www.gin.cnrs.fr/AAL-217?lang=en.

This atlas is the result of an automated anatomical parcellation of the
spatially normalized single-subject high-resolution T1 volume provided by the
Montreal Neurological Institute (MNI) (Collins et al., 1998).

### Creation of atlas

The atlas AAL2 was download from
http://www.gin.cnrs.fr/wp-content/uploads/aal2_for_SPM12.tar.gz on 25.
Sept. 2018. We used the `AAL2.nii` file contained in the `atlas` folder. The
labels are coming from `AAL2.xml` also contain in the `atlas` folder.

### License

Unknown.

## Destrieux 2009

### Background
The 'Destrieux' cortical atlas is based on a parcellation scheme that first
divided the cortex into gyral and sulcal regions, the limit between both being
given by the curvature value of the surface. A gyrus only includes the cortex
visible on the pial view, the hidden cortex (banks of sulci) are marked sulcus.
The result is a complete labeling of cortical sulci and gyri.

### Creation of atlas

The atlas is a direct copy of the file `aparc.a2009s+aseg.mgz` from the folder
`subjects/cvs_avg35_inMNI152/mri` in FreeSurfer version 6.0. The `mgz` file was
converted into NIfTI standard with FreeSurfer's `mri_convert` and the labels
table was created from `FreeSurferColorLUT.txt`.

### License

The Destrieux 2009 atlas is under the License terms of FreeSurfer that can be
found here: https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
**Important Note**: FreeSurfer's license agreement requires you to register. For
more, see this link: https://surfer.nmr.mgh.harvard.edu/fswiki/License

## Desikan & Killiany

The 'Desikan-Killiany' cortical atlas is a gyral based atlas: i.e., a gyrus was
defined as running between the bottoms of two adjacent sulci. That is, a gyrus
includes the part visible on the pial view + adjacent banks of the sulci
limiting this gyrus.

### Creation of atlas

The atlas is a direct copy of the file `aparc+aseg.mgz` from the folder
`subjects/cvs_avg35_inMNI152/mri` in FreeSurfer version 6.0. The `mgz` file was
converted into NIfTI standard with FreeSurfer's `mri_convert` and the labels
table was created from `FreeSurferColorLUT.txt`.

### License

The Destrieux 2009 atlas is under the License terms of FreeSurfer that can be
found here: https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
**Important Note**: FreeSurfer's license agreement requires you to register. For
more, see this link: https://surfer.nmr.mgh.harvard.edu/fswiki/License

## Harvard Oxford Atlas

The Harvard-Oxford cortical and subcortical structural atlases are probabilistic
atlases covering 48 cortical and 21 subcortical structural areas, derived from
structural data and segmentations kindly provided by the
[Harvard Center for Morphometric Analysis](http://www.cma.mgh.harvard.edu/).

T1-weighted images of 21 healthy male and 16 healthy female subjects (ages
18-50) were individually segmented by the CMA using semi-automated tools
developed in-house. The T1-weighted images were affine-registered to MNI152
space using FLIRT (FSL), and the transforms then applied to the individual
labels. Finally, these were combined across subjects to form population
probability maps for each label.

### Creation of atlas

This atlas is a direct copy of the two files
`HarvardOxford-cortl-prob-1mm.nii.gz` and `HarvardOxford-sub-prob-1mm.nii.gz`
Harvard Oxford atlas packaged with FSL 5.0. The files were merged using FSL's
`fslmerge command`.

This probability atlas now also contains the probability values for the labels
`Cerebral_White_Matter` and `Cerebral_Cortex`, spanning the whole brain. Because
the probability for those labels will often be higher than for other labels, we
had to manually remove those labels from the atlas and label file. This can be
done with the following command:

```python
import numpy as np
import nibabel as nb

# Load the image and the data
img = nb.load('atlas_harvard_oxford.nii.gz')
data = img.get_data()

# Remove labels 96, 97, 107 and 108
new_data = np.delete(data, [96, 97, 107, 108], axis=-1)

# Overwrite atlas file
nb.Nifti1Image(new_data, img.affine, img.header).to_filename('atlas_harvard_oxford.nii.gz')
```

### License

The FSL templates and atlases are under the License that can be found
[here](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Licence).

## Jülich histological (cyto- and myelo-architectonic) atlas

A probabilistic atlas created by averaging multi-subject post-mortem cyto- and
myelo-architectonic segmentations, performed by the team of Profs Zilles and
Amunts at the [Research Center Jülich](http://www.fz-juelich.de/inm/inm-1/DE/Home/home_node.html)
and kindly provided by Simon Eickhoff.

The atlas contains 52 grey matter structures and 10 white matter structures.
This is an update to the data used in Eickhoff's
[Anatomy Toolbox](http://www.fz-juelich.de/inm/inm-1/DE/Forschung/_docs/SPMAnatomyToolbox/SPMAnatomyToolbox_node.html) v1.5.
The atlas is based on the miscroscopic and quantitative histological examination
of ten human post-mortem brains. The histological volumes of these brains were
3D reconstructed and spatially normalised into the space of the MNI single
subject template to create a probabilistic map of each area. For the FSL version
of this atlas, these probabilistic maps were then linearly transformed into
MNI152 space.

### Creation of atlas

This atlas is a direct copy of the file `Juelich-prob-1mm.nii.gz` of the Juelich
atlas packaged with FSL 5.0.

### License

The FSL templates and atlases are under the License that can be found
[here](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Licence).


## Talairach atlas

This is a digitised version of the original (coarsely sliced) Talairach atlas
(Lancaster 2000) after the application of a correcting affine transform
(Lancaster 2007) to register it into MNI152 space. For more see http://talairach.org/about.html#Labels.

The atlas was split into two separat atlases. The `talairach_gyurs` atlas
contains the labels listed under **Gyrus** and the `talairach_ba` atlas contains
the labels listed under **Cell type** and the **Brodmann Areas**, as listed
[here](http://www.talairach.org/labels.html).

### Creation of atlas

The talairach atlas was downloaded from http://www.talairach.org/talairach.nii
on the 24. Sept. 2018 using the `fetch_atlas_talairach` function from
`nilearn.datasets`. The two sub-atlases were separated with `nilearn` and the
datatype was converted from `int64` to `int16` with `fslmaths`.

The talairach atlas is identical to the file `Talairach-labels-1mm.nii.gz` that
is packaged with FSL 5.0.

### License

The FSL templates and atlases are under the License that can be found [here](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Licence).

## MarsAtlas (aka MarsAtlas-Colin27-MNI)

The MarsAtlas used in `atlasreader` is the "The MarsAtlas cortical parcellation
of Colin27 in the MNI space". For more information about this atlas see https://meca-brain.org/software/marsatlas-colin27/.

### Creation of atlas

The atlas was downloaded from https://meca-brain.org/software/marsatlas-colin27/ on the 15.
Oct. 2018 via the following link: https://www.dropbox.com/s/ndz8qtqblkciole/MarsAtlas-MNI-Colin27.zip?dl=1

The labels file was manually created according to the table of cortical(http://meca-brain.org/software/marsatlas/)
and sub-cortical (http://meca-brain.org/software/marsatlas-subcortical/) regions.
For cortical regions, the label name is a combination of the "Full Name" and the
"Brodman Area" label.

### Licesne

For more information about the license and copyright of the `MarsAtlas-Colin27-MNI`
atlas, go to https://meca-brain.org/software/marsatlas-colin27/.


## Neuromorphometrics

Maximum probability tissue labels derived from the "MICCAI 2012 Grand Challenge
and Workshop on Multi-Atlas Labeling" (https://masi.vuse.vanderbilt.edu/workshop2012/index.php/Challenge_Details).
These data were released under the Creative Commons Attribution-NonCommercial
(CC BY-NC) with no end date. Users should credit the MRI scans as originating
from the OASIS project (http://www.oasis-brains.org/) and the labeled data as
"provided by Neuromorphometrics, Inc. (http://Neuromorphometrics.com/) under
academic subscription". These references should be included in all workshop and
final publications.

### Creation of atlas

The Neuromorphometrics atlas is a direct copy of the file `labels_Neuromorphometrics.nii`
in the folder `spm12/tpm`, available within SPM12.

### Licesne

The direct license of the atlas can be found here: http://www.neuromorphometrics.com/wp-content/uploads/2013/06/NVM_Demo_License2013_v1.txt

The Neuromorphometrics atlas version that we are using is the one provided with
SPM12, therefore we also want to point to the SPM license of SPM12.
As SPM is free but copyright software, distributed under the terms of the GNU
General Public Licence as published by the Free Software Foundation (either
version 2, as given in file spm_LICENCE.man, or at your option, any later
version). Further details on "copyleft" can be found at http://www.gnu.org/copyleft/.
In particular, SPM is supplied as is. No formal support or maintenance is
provided or implied. For more see: https://www.fil.ion.ucl.ac.uk/spm/software/
or https://github.com/neurodebian/spm12/blob/master/LICENCE.txt

## Atlas of Intrinsic Connectivity of Homotopic Areas (AICHA)

AICHA (Atlas of Intrinsic Connectivity of Homotopic Areas) is a functional brain
ROIs atlas based on resting-state fMRI data acquired in 281 individuals.
AICHA ROIs cover the whole cerebrum, each having 1. homogeneity of its
constituting voxels intrinsic activity, and 2. a unique homotopic contralateral
counterpart with which it has maximal intrinsic connectivity.

The AICHA atlas includes 192 couples of homotopic regions for a total of 384
regions. AICHA is provided in the MNI stereotaxic space (MNI ICBM 152, Template
sampling size of 2x2x2 mm3 voxels; bounding box, x = -90 to 90 mm, y = -126 to
91 mm, z = -72 to 109 mm). Each region get a pseudo-color with odd number for
region belong to the left hemisphere and even for the right. Each homotopic pair
is labeled with an odd number (Left) and the following even number (Right).
For example: 1 and 2 code for G_Frontal_Sup-1-L and G_Frontal_Sup-1-R
respectively, 3 and 4 code for G_Frontal_Sup-2-L and G_Frontal_Sup-2-R,...

AICHA atlas includes both regions located in the crown of the gyri (named Gyrus,
region name beginning by "G_") and regions located in the depth of the sulci
(named Suclus, region name beginning by "S_"). In addition the subcortical
nuclei were labeled separately (name Nucleus, region name beginning by "N_").
Different parcels belonging to the same anatomical region were labeled with
numbers (starting to 1). For example the precuneus show as 9 subparts labeled
from G_Precuneus-1-L to G_Precuneus-9-L.

### Creation of atlas

The atlas was downloaded from http://www.gin.cnrs.fr/en/tools/aicha/ on the 24.
Sept. 2018 via the following link: http://www.gin.cnrs.fr/wp-content/uploads/aicha_v1.zip

We compressed the file `AICHA.nii` with `fslmaths` and transformed the data type
from `float32` to `int16`.

### License

This atlas is protected by copyright; you can freely use it for none profit
research purposes, providing the above reference is cited. For other use please
contact us through aicha.gin.brainatlas@gmail.com

# References

If you're using `atlasreader`, please make sure to cite the references that
correspond to the atlases you are using:



## MNI152 T1 1mm Template

The references are related to FSL 5.0:

- Rahul S. Desikan, Florent Ségonne, Bruce Fischl, Brian T. Quinn, Bradford C. Dickerson, Deborah Blacker, Randy L. Buckner, Anders M. Dale, R. Paul Maguire, Bradley T. Hyman, Marilyn S. Albert, Ronald J. Killiany (2006). An automated labeling system for subdividing the human cerebral cortex on MRI scans into gyral based regions of interest. Neuroimage, 31, 968–980.
- Simon B Eickhoff, Klaas E Stephan, Hartmut Mohlberg, Christian Grefkes, Gereon R Fink, Katrin Amunts, Karl Zilles (2005). A new SPM toolbox for combining probabilistic cytoarchitectonic maps and functional imaging data. Neuroimage, 25, 1325–1335.
- Kegang Hua, Jiangyang Zhang, Setsu Wakana, Hangyi Jiang, Xin Li, Daniel S Reich, Peter A Calabresi, James J Pekar, Peter C M van Zijl, Susumu Mori (2008). Tract probability maps in stereotaxic spaces: analyses of white matter anatomy and tract-specific quantification. Neuroimage, 39, 336–347.
- J Mazziotta, A Toga, A Evans, P Fox, J Lancaster, K Zilles, R Woods, T Paus, G Simpson, B Pike, C Holmes, L Collins, P Thompson, D MacDonald, M Iacoboni, T Schormann, K Amunts, N Palomero-Gallagher, S Geyer, L Parsons, K Narr, N Kabani, G Le Goualher, D Boomsma, T Cannon, R Kawashima, and B Mazoyer (2001). A probabilistic atlas and reference system for the human brain: International Consortium for Brain Mapping (ICBM). The Royal Society Philosophical Transactions B, 356, 1293–1322.
- Jack L. Lancaster, Marty G. Woldorff, Lawrence M. Parsons, Mario Liotti, Catarina S. Freitas, Lacy Rainey, Peter V. Kochunov, Dan Nickerson, Shawn A. Mikiten, Peter T. Fox (2000). Automated Talairach atlas labels for functional brain mapping. Human Brain Mapping, 10, 120–131.
- TEJ Behrens, H. Johansen-Berg, MW Woolrich, SM Smith, CAM Wheeler-Kingshott, PA Boulby, GJ Barker, EL Sillery, K. Sheehan, O. Ciccarelli , AJ Thompson, JM Brady, PM Matthews (2003). Non-invasive mapping of connections between human thalamus and cortex using diffusion imaging. Nature Neuroscience, 6, 750–757.
- Jörn Diedrichsen, Joshua H Balsters, Jonathan Flavell, Emma Cussans, Narender Ramnani (2009). A probabilistic MR atlas of the human cerebellum. Neuroimage, 46, 39-46.

## ICBM 2009c Nonlinear Asymmetric

- VS Fonov, AC Evans, K Botteron, CR Almli, RC McKinstry, DL Collins and BDCG, Unbiased average age-appropriate atlases for pediatric studies, NeuroImage,Volume 54, Issue 1, January 2011, ISSN 1053–8119, DOI: 10.1016/j.neuroimage.2010.07.033
- VS Fonov, AC Evans, RC McKinstry, CR Almli and DL Collins, Unbiased nonlinear average age-appropriate brain templates from birth to adulthood, NeuroImage, Volume 47, Supplement 1, July 2009, Page S102 Organization for Human Brain Mapping 2009 Annual Meeting, DOI: http://dx.doi.org/10.1016/S1053-8119(09)70884-5

## Anatomical Automatic Labeling 2 (AAL2)

- Automated Anatomical Labeling of Activations in SPM Using a Macroscopic Anatomical Parcellation of the MNI MRI Single-Subject Brain. N. Tzourio-Mazoyer, B. Landeau, D. Papathanassiou, F. Crivello, O. Étard, N. Delcroix, B. Mazoyer, and M. Joliot. NeuroImage 2002. 15 :273-289
http://dx.doi.org/10.1006/nimg.2001.0978
- Implementation of a new parcellation of the orbitofrontal cortex in the automated anatomical labeling atlas. Rolls ET, Joliot M & Tzourio-Mazoyer N (2015) . NeuroImage
http://dx.doi.org/10.1016/j.neuroimage.2015.07.075
- Collins, D. L., Zijdenbos, A. P., Kollokian, V., Sled, J. G., Kabani, N. J., Holmes, C. J., & Evans, A. C. (1998). Design and construction of a realistic digital brain phantom. IEEE transactions on medical imaging, 17(3), 463-468.

## Destrieux 2009

- Fischl, Bruce, et al. "Automatically parcellating the human cerebral cortex." Cerebral cortex 14.1 (2004): 11-22.
- Destrieux, C., et al. "A sulcal depth-based anatomical parcellation of the cerebral cortex." NeuroImage 47 (2009): S151.

## Desikan & Killiany

- Fischl, Bruce, et al. "Automatically parcellating the human cerebral cortex." Cerebral cortex 14.1 (2004): 11-22.
- An automated labeling system for subdividing the human cerebral cortex on MRI scans into gyral based regions of interest, Desikan et al., (2006). NeuroImage, 31(3):968-80.

## Harvard Oxford Atlas

- Makris N, Goldstein JM, Kennedy D, Hodge SM, Caviness VS, Faraone SV, Tsuang MT, Seidman LJ. Decreased volume of left and total anterior insular lobule in schizophrenia. Schizophr Res. 2006 Apr;83(2-3):155-71
- Frazier JA, Chiu S, Breeze JL, Makris N, Lange N, Kennedy DN, Herbert MR, Bent EK, Koneru VK, Dieterich ME, Hodge SM, Rauch SL, Grant PE, Cohen BM, Seidman LJ, Caviness VS, Biederman J. Structural brain magnetic resonance imaging of limbic and thalamic volumes in pediatric bipolar disorder. Am J Psychiatry. 2005 Jul;162(7):1256-65
- Desikan RS, Ségonne F, Fischl B, Quinn BT, Dickerson BC, Blacker D, Buckner RL, Dale AM, Maguire RP, Hyman BT, Albert MS, Killiany RJ. An automated labeling system for subdividing the human cerebral cortex on MRI scans into gyral based regions of interest. Neuroimage. 2006 Jul 1;31(3):968-80.
- Goldstein JM, Seidman LJ, Makris N, Ahern T, O'Brien LM, Caviness VS Jr, Kennedy DN, Faraone SV, Tsuang MT. Hypothalamic abnormalities in schizophrenia: sex effects and genetic vulnerability. Biol Psychiatry. 2007 Apr 15;61(8):935-45

## Jülich histological (cyto- and myelo-architectonic) atlas

- Eickhoff et al., A new SPM toolbox for combining probabilistic cytoarchitectonic maps and functional imaging data. Neuroimage 25(4):1325-35 (2005)
- Eickhoff et al., Testing anatomically specified hypotheses in functional imaging using cytoarchitectonic maps. NeuroImage 32(2): 570-582 (2006)
- Eickhoff et al., Assignment of functional activations to probabilistic cytoarchitectonic areas revisited. NeuroImage, 36(3): 511-521 (2007)

## Talairach atlas

- Talairach et al. Co-planar stereotaxic atlas of the human brain. Thieme, New York. (1988)
- Lancaster et al. Bias between MNI and Talairach coordinates analyzed using the ICBM-152 brain template. Human Brain Mapping (in press) (2007)
- Lancaster JL, Woldorff MG, Parsons LM, Liotti M, Freitas CS, Rainey L, Kochunov PV, Nickerson D, Mikiten SA, Fox PT, “Automated Talairach Atlas labels for functional brain mapping”. Human Brain Mapping 10:120-131, 2000.
- Lancaster JL, Rainey LH, Summerlin JL, Freitas CS, Fox PT, Evans AC, Toga AW, Mazziotta JC. Automated labeling of the human brain: A preliminary report on the development and evaluation of a forward-transform method. Hum Brain Mapp 5, 238-242, 1997.

## MarsAtlas

For information about references, see https://meca-brain.org/software/marsatlas-colin27/.

## Neuromorphometrics

Reference missing.

## Atlas of Intrinsic Connectivity of Homotopic Areas (AICHA)

- Joliot M, Jobard G, Naveau M, Delcroix N, Petit L, Zago L, Crivello F, Mellet E, Mazoyer B, Tzourio-Mazoyer N (2015) AICHA: An atlas of intrinsic connectivity of homotopic areas. J Neurosci Methods 254:46-59.
