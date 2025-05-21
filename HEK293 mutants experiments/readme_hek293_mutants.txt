Centromere-Tau Colocalization and Clustering Analysis
Last Updated: 2025.05.21

This repository contains MATLAB scripts for analyzing colocalization between centromere, tau, and tubulin signals in fixed-cell confocal images. The analysis includes quantitative image extraction, scatter-based clustering, and comparative analysis across tau variants.

-------------------------
1. Script Overview
-------------------------

- centromere_analysis.m (Prometaphase)
  - Loads multi-channel confocal image stacks.
  - Extracts and aligns centromere-centered regions per cell.
  - Computes Pearson correlation coefficients (PCC) between tau, tubulin, and centromere channels (raw, cropped, and intensity profiles).
  - Visualizes representative cells, intensity profiles, and cross-channel scatter plots.
  - Produces detailed panel figures and summary statistics for bar plots and publication figures.

- analysis2.m (Prometaphase)
  - Compiles and compares Pearson correlation coefficients (PCCs) across wild-type and mutant tau conditions.
  - Creates boxplots and spread plots for:
  - Tubulin-Tau correlation
  - Centromere-Tau correlation
  - Group labels: Wild Type, T231D/S235D, S262D.

- clustering_line.m (Metaphase)
  - Performs tau-DNA scatter plot analysis in a selected ROI.
  - Allows manual definition of decision boundaries (via `drawline`) to define a custom coordinate frame.
  - Applies cosine-distance k-means clustering to tau-DNA data.
  - Identifies cluster with maximum slope (indicative of strong co-alignment).
  - Generates filtered cross-correlation maps and joint histograms.

-------------------------
2. Input Requirements
-------------------------
- Multi-channel TIFF files (channels: Tubulin, Tau, Centromere, DNA).
- Directory structure organized by condition (e.g., `WT\`, `Double\`, `S262D\`).
- Optional: ROI selection and registration metadata (automatically generated).

-------------------------
3. Output Files
-------------------------
- `img_crop`, `img_centromere`, `I_centromere`: image crops, aligned stacks, and intensity profiles.
- `PCC`, `nmol`: correlation coefficients and molecule count per image.
- `.fig` plots and `.mat` files for all results.
- Final summary plots: scatter plots, kymographs, bar graphs, cross-correlation maps, and RGB channel overlays.

-------------------------
4. Dependencies
-------------------------
- MATLAB (R2021a or newer)
- Signal Processing Toolbox
- Custom or helper functions:
- loadTifStack16, imshow2, maxfig, subplot2
- countSM (molecule detection), centroid (ROI alignment)
- plotSpread (scatter+summary), rescale, msum, mmean
- clustering_line: requires user input via `drawline`

-------------------------
5. Contact
-------------------------
For questions or contributions, please contact:  
📧 mjshon@postech.ac.kr
