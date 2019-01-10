"""
This script crops the atlases in the atlases folder, changes their datatype to
a minimum and compresses the files with a high enough compression level. This
all to reduce the disk size and control the load time (driven by the level of
compression).
"""

import nibabel as nb
from glob import glob
from nilearn.image import crop_img

# Get the list of atlases
atlases = sorted(glob('atlases/atlas_*nii.gz'))

for fname in atlases:

    # Load atlas and crop image
    img = crop_img(fname)

    # Get data array
    data = img.get_data()

    # Decide which datatype to use
    if data.max() <= 255 and data.min() >= 0:
        dtype = '>u1'
    elif data.max() <= 65535 and data.min() >= 0:
        dtype = '>u2'
    else:
        dtype = 'i2'

    # Create a new image with the correct datatype
    img.set_data_dtype(dtype)
    new_img = nb.Nifti1Image(
        data.astype(dtype), affine=img.affine, header=img.header)

    # Change the compression level of the NIfTI image
    nb.openers.Opener.default_compresslevel = 6

    # Overwrite previous atlas file with new image
    new_fname = fname.replace('.nii.gz', '.nii.gz')
    new_img.to_filename(new_fname)
