Harvard-Oxford atlas
--------------------

The Harvard-Oxford cortical and subcortical structural atlases are a copy of
FSL's (Version 5.0) under `$FSL_DIR/data/atlases/HarvardOxford/`:
 - HarvardOxford-cortl-prob-1mm.nii.gz' 
 - HarvardOxford-sub-prob-1mm.nii.gz' 
For further information see http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Atlases.

The following 4 areas were deleted from the subcortical structural atlas, as
they are already coverd in the cortical structural atlas:
 - 'Left Cerebral White Matter'
 - 'Left Cerebral Cortex '
 - 'Right Cerebral White Matter'
 - 'Right Cerebral Cortex '

Afterwards, the cortical and subcortical structural atlases were merged with
the following python code:

    import nibabel as nb
    cort = nb.load('HarvardOxford-cort-prob-1mm.nii.gz')
    sub = nb.load('HarvardOxford-sub-prob-1mm.nii.gz')
    newsubdata = sub.get_data()[:, :, :, range(2, 11) + range(13, 21)]
    data = np.concatenate((cort.get_data(), newsubdata), axis=3)
    affine = cort.get_affine()
    nb.Nifti1Image(data.astype('uint8'), affine).to_filename('bladi.nii.gz')
 

References
----------

Makris N, Goldstein JM, Kennedy D, Hodge SM, Caviness VS, Faraone SV, Tsuang MT, Seidman LJ. Decreased volume of left and total anterior insular lobule in schizophrenia. Schizophr Res. 2006 Apr;83(2-3):155-71

Frazier JA, Chiu S, Breeze JL, Makris N, Lange N, Kennedy DN, Herbert MR, Bent EK, Koneru VK, Dieterich ME, Hodge SM, Rauch SL, Grant PE, Cohen BM, Seidman LJ, Caviness VS, Biederman J. Structural brain magnetic resonance imaging of limbic and thalamic volumes in pediatric bipolar disorder. Am J Psychiatry. 2005 Jul;162(7):1256-65

Desikan RS, Ségonne F, Fischl B, Quinn BT, Dickerson BC, Blacker D, Buckner RL, Dale AM, Maguire RP, Hyman BT, Albert MS, Killiany RJ. An automated labeling system for subdividing the human cerebral cortex on MRI scans into gyral based regions of interest. Neuroimage. 2006 Jul 1;31(3):968-80.

Goldstein JM, Seidman LJ, Makris N, Ahern T, O'Brien LM, Caviness VS Jr, Kennedy DN, Faraone SV, Tsuang MT. Hypothalamic abnormalities in schizophrenia: sex effects and genetic vulnerability. Biol Psychiatry. 2007 Apr 15;61(8):935-45
