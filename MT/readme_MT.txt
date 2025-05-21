Magnetic Tweezers Data Analysis for WLC and Force Clamp Experiments

Last Updated: 2025.05.21

This repository contains MATLAB scripts for analyzing single-molecule force spectroscopy experiments using magnetic tweezers. 
The analysis is divided into two modes: Worm-Like Chain (WLC) stretching experiments and Force Clamp (FC) measurements. 
Each mode has a data preprocessing and a main analysis script.

===================================
1. WLC Analysis Scripts
===================================

- analysis_WLC.m
  - Preprocesses raw magnetic tweezers data from WLC experiments.
  - Loads tracking, force calibration, and positional data for multiple beads.
  - Computes corrected bead displacements and interpolated force/piezo traces.
  - Saves intermediate data to `analysis_WLC.mat`.

- analysis2_WLC.m
  - Performs full analysis and visualization of WLC data.
  - Identifies force ranges for stretching and condensation phases.
  - Compares experimental force-extension curves (FECs) to theoretical WLC model.
  - Extracts unzipping and re-zipping forces across conditions (e.g., tau concentration).
  - Produces bar plots summarizing results.

===================================
2. Force Clamp (FC) Analysis Scripts
===================================

- analysis_FC.m
  - Preprocesses raw data from force clamp experiments.
  - Similar structure to WLC preprocessing; computes corrected coordinates and force.
  - Saves processed data into `analysis_FC.mat`.

- analysis2_FC.m
  - Analyzes displacement trajectories under constant force conditions.
  - Identifies event timings and calculates displacement derivatives (velocity).
  - Extracts step sizes using peak detection in smoothed derivative signals.
  - Compiles step distributions and generates histograms.

===================================
Data Requirements
===================================
Each experiment requires:
- `rXXX-XXX.xls` : Bead position data (magnetic beads + reference beads).
- `sXXX-XXX.xls` : Magnet position/force metadata.
- `cXXX.fps`     : Calibration info (sampling rate, offsets, orientation).

Raw data should be organized in:
- `WLC/raw*/` for WLC experiments.
- `Force clamp/raw*/` for FC experiments.

===================================
Dependencies
===================================
- MATLAB with Signal Processing Toolbox
- Custom functions (e.g., `WLC_inv`, `correctFEC`, `closest`, `vline`, `maxfig`)
  [Ensure these functions are included or available in your path]


===================================
Contact
===================================
For questions or contributions, please contact: mjshon@postech.ac.kr
