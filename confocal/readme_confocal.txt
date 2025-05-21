Centromere Imaging Analysis for Cell Cycle Stages
Last Updated: 2025.05.21

This folder contains a MATLAB script for analyzing 4-channel confocal images of HEK293 and SH-SY5Y cells during mitosis. 
The focus is on centromere signal extraction and visualization at prometaphase, metaphase, and anaphase.

-------------------------
1. Main Script
-------------------------
- analysis.m
  - Loads TIFF image stacks across experimental folders.
  - Extracts centromere-centered regions using coordinates (cnt).
  - Normalizes fluorescence channels and overlays RGB images.
  - Aligns centromeres and extracts intensity profiles along the x-axis.
  - Computes Pearson correlation (PCC) between selected channels.
  - Generates figures: full-cell views, aligned centromeres, profile plots.

-------------------------
2. Required Data
-------------------------
- .tif image stacks (4-channel) under experiment folders.
- cnt variable: [cell_index, x, y] for each centromere (manual or preloaded).

-------------------------
3. Output Files
-------------------------
- analysis2.mat: cropped images, aligned centromeres, PCC values, intensity profiles.
- images.fig: full-cell views
- centromere.fig: aligned centromere panels
- results.fig: summary plots (mean images, foci count, PCC, intensity profiles)

-------------------------
4. Dependencies
-------------------------
- MATLAB (R2021 or newer)
- Custom functions required in path:
  loadTifStack16, centroid, maxfig, sfigure, subplot2, imshow2, plotSpread, etc.

-------------------------
5. Contact
-------------------------
For questions, contact: mjshon@postech.ac.kr
