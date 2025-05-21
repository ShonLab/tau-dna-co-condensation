Low-Magnification Imaging Analysis of Tau Effects
Last Updated: 2025.05.21

This script analyzes 2-channel low-magnification confocal images comparing conditions with and without tau. It includes image registration, crosstalk correction, and line-like structure detection.

-------------------------
1. Main Script
-------------------------
- analysis_low mag.m
  - Loads and processes TIFF image stacks (2 channels per field).
  - Applies affine image registration using manually selected landmarks.
  - Corrects channel crosstalk via robust linear regression.
  - Extracts cropped regions of interest (ROI) for visualization.
  - Detects and counts linear structures (e.g., MT bundles) using skeletonization and filtering.
  - Merges nearby segments with similar orientation.
  - Generates summary plots and saves results.

-------------------------
2. Input Requirements
-------------------------
- TIFF images in `raw2_low mag\` and `raw3_low mag\` directories.
  - Channel 1: tubulin (or similar)
  - Channel 2: DNA or tau label
- Manual control point registration (via cpselect) stored in `registration.mat`.

-------------------------
3. Output Files
-------------------------
- analysis_raw2.mat / analysis_raw3.mat : loaded and preprocessed image data
- analysis2.mat : processed image crops, segment counts (`nline`)
- Figures: visualization of raw/corrected channels and detected line structures

-------------------------
4. Dependencies
-------------------------
- MATLAB (R2021 or newer)
- Signal Processing Toolbox
- Custom functions required in path:
  loadTifStack16, imshow2, maxfig, subplot2, rescale

-------------------------
5. Contact
-------------------------
For questions, contact: mjshon@postech.ac.kr
