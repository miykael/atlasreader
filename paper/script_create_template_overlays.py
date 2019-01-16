from nilearn import plotting, image
from glob import glob
import nibabel as nb
import numpy as np

atlases = glob('atlases/*gz')
template = 'templates/mni_icbm152_t1_tal_nlin_asym_09c_brain.nii.gz'

for a in atlases:
    atlas_dim = len(nb.load(a).shape)

    if atlas_dim == 4:
        plotting.plot_prob_atlas(
            a, bg_img=template, title=a[15:-7],  cut_coords=[10,0,0],
            threshold=0.5, draw_cross=False, output_file=a[15:-7]+'.png')
    else:

        # Shuffle color values to get better disceprancy between neighbors
        img = image.load_img(a)
        data = img.get_data()
        orig = np.unique(data)[1:]
        rand = np.arange(1, len(orig) + 1)
        np.random.shuffle(rand)
        for i, o in enumerate(orig):
            data[data==o] = rand[i]
        data[data!=0] -= 10000
        img = image.new_img_like(img, data)

        plotting.plot_roi(
            img, bg_img=template, title=a[15:-7], cut_coords=[10,0,0],
            draw_cross=False, output_file=a[15:-7]+'.png')
