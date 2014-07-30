fmri_ica_artifactremoval
========================

A Matlab toolbox for fMRI artifact reduction using ICA and global decision trees.

A matlab package for the artifact identification using Global Decision Tree classifiers and for their training.
To use this package you will also need NiFTI tools 

http://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image

and FSL software package

http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/

(Or, more precisely, Melodic 3.0 is needed). 

This toolbox is designed to run under Linux or Unix.

Installation by just unpacking the contents of the zip file. Running the functions might 
require that fsl binaries are available on the PATH. Manuals subdir contains manuals that refer 
to an older version of the toolbox (1.1 for fmri_ica_classify (current 1.1.2) and 1.2 for fmri_ica_training (current 1.2.1)). 
In particular, Melodic version 3.0 is now expected (instead of Melodic 2.0). 

Readily trained classifiers are provided in classfiers subdir (please read the manual) and example training 
file is on example_training_data subdir.

Reference:

J. Tohka , K. Foerde, A.R. Aron, S. M. Tom, A.W. Toga, and R.A. Poldrack. 
Automatic Independent Component Labeling for Artifact Removal in fMRI. 
NeuroImage, 39(3):1227-45, 2008.
http://www.ncbi.nlm.nih.gov/pubmed/18042495


