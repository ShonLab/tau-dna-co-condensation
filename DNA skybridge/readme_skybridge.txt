Single-Molecule Colocalization & Diffusion Analysis (Tau-DNA Imaging)
Last Updated: 2025.05.21

This repository contains MATLAB scripts for analyzing dual-channel time-lapse fluorescence imaging data, focusing on tau and DNA colocalization. The workflow includes background subtraction, vertical artifact masking, affine registration, particle detection, tracking, diffusion coefficient calculation, and result visualization (including kymographs and bar plots).

-------------------------
1. Main Script
-------------------------

- analysis.m
  - Loads dual-channel TIFF images for multiple experimental conditions.
  - Subtracts background and detects vertical line artifacts (e.g., scanline noise).
  - Applies affine registration to align tau and DNA channels.
  - Identifies colocalized spots and tracks them over time.
  - Calculates diffusion coefficients (D) and anomalous exponents (α).
  - Computes tau-DNA correlation coefficients using kymographs.
  - Saves composite images, trajectories, and summary figures.
  - Plots bar graphs for colocalization ratio, log(D), α, and CCF.

-------------------------
2. Input Requirements
-------------------------
- Dual-channel .tif images (tau = top half, DNA = bottom half).
- Folder structure: subfolders per condition, excluding any named "result".
- A manually defined affine transform saved as `tform.mat` (optional; created on first run).

-------------------------
3. Output Files
-------------------------
- result/: composite images, kymographs, trajectory overlays
- data.mat: stores all tracking and diffusion data (used in plotting section)
- Bar plots and histograms:
  - Colocalization ratio
  - Log(D), α (alpha), and correlation coefficient (CCF) for colocalized vs. non-colocalized spots

-------------------------
4. Dependencies
-------------------------
- MATLAB (R2021a or newer)
- Signal Processing Toolbox
- Required custom functions:
  - `vertical_line_detection.m`: masks periodic vertical line artifacts
  - `dualviewer_merger.m`: splits and registers dual-view images
  - `colocal_1d.m`: detects colocalized spots between two channels
  - `countSM.m`: detects single molecules via local maxima and centroid filtering
  - `particle_tracking.m`: tracks particles frame-to-frame using nearest neighbor logic
  - `CalD_1D.m`: computes diffusion coefficient (D) and anomalous exponent (α) via MSD analysis
  - `kymo.m`: extracts kymograph slices over time
  - `trajectory.m`: generates 2D binary maps of particle trajectories
  - `bresenham.m`: efficient line rasterization between two points

-------------------------
5. Contact
-------------------------
For questions, suggestions, or contributions, please contact:
📧 mjshon@postech.ac.kr
