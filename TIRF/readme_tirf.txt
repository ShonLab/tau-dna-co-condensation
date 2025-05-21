DNA Tethering Experiments - Image Analysis Scripts
Last Updated: 2025.05.21

This repository contains MATLAB scripts for analyzing fluorescence imaging of tau-DNA interactions using single and doubly tethered DNA constructs. The analysis covers signal extraction, transformation, alignment, dynamic quantification, and visualization.

-------------------------
1. Scripts Overview
-------------------------

- analysis_doubly tethered DNA.m
  - Processes time-lapse dual-channel TIFFs from multi-condition experiments.
  - Computes background-corrected averages, min/max projections.
  - Registers tau channel to DNA using affine transforms.
  - Identifies tether endpoints and detects connected molecules.
  - Aligns molecules along their axis and extracts fluorescence profiles.
  - Measures envelope width (via sinusoidal fitting) and correlation with tau.
  - Outputs include: averaged images, ROI movies, fluorescence plots, kymographs, and envelope fits.

- analysis_doubly tethered DNA_continued.m
  - Compiles summary figures across all conditions.
  - Plots envelope width vs end-to-end distance and tau concentration.
  - Generates bar plots for intensity, roughness, and Pearson correlations.
  - Produces merged movies for representative molecules.

- analysis_singly tethered DNA.m
  - Handles registration and cropping for single tethered molecules.
  - Extracts RGB overlay frames and compiles AVI movie.
  - Generates kymographs and estimates DNA length changes over time.

-------------------------
2. Input Requirements
-------------------------
- Dual-channel TIFF files (one for DNA, one for tau).
- Directory structure: `raw*/` per condition.
- Registration control points stored via `cpselect` (or interactive).
- Predefined masks and thresholds are hardcoded; adjust as needed.

-------------------------
3. Output Files
-------------------------
- analysis.mat / analysis2.mat : processed image data and metadata.
- *.jpg / *.avi : raw and processed images, montages, movies, and summary plots.
- Kymographs and envelope width results for each molecule.

-------------------------
4. Dependencies
-------------------------
- MATLAB (R2021 or newer)
- Signal Processing Toolbox
- Custom functions required in path:
  loadTifStack16, imshow2, maxfig, plotSpread, mmean, countSM, rescale, linePixels, cpselect (GUI-based), fitgeotform2d

-------------------------
5. Contact
-------------------------
For questions, contact: mjshon@postech.ac.kr
