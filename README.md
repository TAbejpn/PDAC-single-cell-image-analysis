## Files
- `Ecad_GATA6_S100A2_Analysis_public.m` - main MATLAB analysis script
- `checkbbox_v2.m` - helper function for bounding-box adjustment
- `example_HE.jpg` - representative H&E image
- `example_FL.png` - three-channel fluorescence image
- `example_GATA6.tif` - GATA6 fluorescence image
- `example_HoverNet.mat` - HoVer-Net nuclear segmentation

The fluorescence image contains the following channels:
1. S100A2
2. E-cadherin (CDH1)
3. DAPI

## Trained model
The trained ResNet50 model (`example_Ecad_ResNet50.mat`) required for CDH1 pattern classification is available from Zenodo:

https://doi.org/10.5281/zenodo.24387393

Version: 1.0.0

## Requirements
The analysis was implemented in MATLAB.
Required MATLAB functionality includes the Image Processing Toolbox, Computer Vision Toolbox, and Deep Learning Toolbox. ResNet50 support is also required.

## Usage
Download `example_Ecad_ResNet50.mat` from Zenodo and place it in the same directory as the MATLAB scripts and example files.
Run:
`Ecad_GATA6_S100A2_Analysis_public.m`
The script performs cell-level classification and saves annotated images in the `output` subfolder.


## Archived version
The released version of this repository is archived in Zenodo:
https://doi.org/10.5281/zenodo.24389982
