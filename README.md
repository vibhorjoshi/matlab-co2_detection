# 🌍 CO₂ Detection from Hyperspectral Imagery using AVIRIS Data

<p align="center">
  <img src="https://img.shields.io/badge/MATLAB-R2023b%2B-blue.svg" alt="MATLAB Version">
  <img src="https://img.shields.io/badge/Challenge-MathWorks_Project_251-orange.svg" alt="MathWorks Challenge">
  <img src="https://img.shields.io/badge/License-MIT-green.svg" alt="License">
</p>

## 📌 Project Overview

This repository presents a comprehensive, physics-aware MATLAB framework for the detection, visualization, and statistical analysis of atmospheric CO₂ signatures using **Airborne Visible/Infrared Imaging Spectrometer (AVIRIS)** hyperspectral imagery.

Developed as a formal submission for the **MATLAB and Simulink Challenge (Project 251)**, this project directly addresses the physical and mathematical limitations of standard anomaly detectors. Standard Matched Filters routinely fail when applied to raw reflectance data because gaseous anomalies present as absorption features (negative dips) rather than emissive features (positive spikes). 

To solve this, this repository implements a **Progressive Spectral Conditioning Pipeline** that evaluates baseline continuum removal techniques against a newly proposed **Cluster-Tuned Adaptive Coherence Estimator (CT-ACE)**. By combining the Beer-Lambert physical absorbance transformation with scale-invariant localized covariance mathematics, this framework achieves near-perfect background clutter suppression and zero false alarms.

---

## 📑 Table of Contents
1. [Theoretical Background: The Bottleneck](#1-theoretical-background-the-bottleneck)
2. [The Progressive Pipeline Architecture](#2-the-progressive-pipeline-architecture)
3. [Visual and Quantitative Results](#3-visual-and-quantitative-results)
   - [Stagewise Ablation & Spatial Isolation](#stagewise-ablation--spatial-isolation)
   - [Profile Analysis & Noise Floor](#profile-analysis--noise-floor)
   - [Threshold Sensitivity & False Alarms](#threshold-sensitivity--false-alarms)
   - [Selectivity Metrics (SBR)](#selectivity-metrics-sbr)
   - [Geospatial Validation](#geospatial-validation)
4. [Toolbox Utilization](#4-toolbox-utilization)
5. [Data Acquisition & Setup](#5-data-acquisition--setup)
6. [Repository Structure & Usage](#6-repository-structure--usage)
7. [Citation & License](#7-citation--license)

---

## 1. Theoretical Background: The Bottleneck

### The Physical Flaw in Standard Filters
Most hyperspectral target detection algorithms (such as the Spectral Angle Mapper or Global Matched Filter) were designed to find solid objects (e.g., vehicles, minerals) that reflect light. Atmospheric gases like CO₂, however, *absorb* specific wavelengths of shortwave infrared (SWIR) light. 

When standard algorithms search for a gas template within a raw reflectance datacube, the mathematical dot-product projects a negative score. To display an image, researchers typically clip negative values to zero, which inadvertently erases the actual gas plume. 

### The Mathematical Flaw (Scale Variance)
Furthermore, a standard global matched filter calculates a single covariance matrix for the entire landscape. This causes the filter to be "scale-variant." If a scene contains deep shadows, highly reflective concrete, and dark vegetation, the standard filter cannot adapt to the local variance, triggering massive false-positive anomalies over bright surfaces.

---

## 2. The Progressive Pipeline Architecture

To systematically prove the necessity of our proposed solution, the framework evaluates four distinct algorithms, treating CO₂ retrieval as a progressive spectral conditioning operation.

### Baseline 1: Continuum Interpolated Band Ratio (CIBR)
CIBR serves as a rudimentary, computationally inexpensive anomaly detector. It measures the relative depth of the known CO₂ absorption feature near **2.05 μm** by interpolating a baseline between a left continuum (2000–2020 nm) and a right continuum (2080–2100 nm). While it highlights the general area of the gas, it is severely compromised by surface albedo variations, resulting in high background noise.

### Baseline 2: Joint Reflectance and Gas Estimator (JRGE)
JRGE attempts to mathematically remove the background continuum using spline-based estimation. By subtracting the continuum slope, it suppresses broadband reflectance variations and mitigates horizontal striping artifacts. However, it still struggles to differentiate the weak gas signal from complex, heterogeneous terrestrial surfaces.

### Baseline 3: Spectral Fitting Algorithm (SFA)
SFA expands the analysis to the entire **1500–2100 nm SWIR region**. It utilizes a synthesized dual-Gaussian CO₂ template centered at 1575 nm and 2005 nm. It projects the raw spectral data onto this template. Because it does not account for the specific covariance of the background interference, its Signal-to-Background ratio remains poor.

### Proposed Solution: Cluster-Tuned Adaptive Coherence Estimator (CT-ACE)
The proposed CT-ACE framework completely solves the limitations of the baselines through a rigorous physics-and-math approach:

1. **Beer-Lambert Absorbance Transform:** Before any filtering occurs, the raw Reflectance ($R$) is converted to Absorbance ($A$) using the relationship:
   $$A = -\log(R)$$
   In absorbance space, gaseous absorption features become positive, additive spikes. This immediately resolves the negative-clipping bug.
2. **Water Vapor Masking:** The opaque atmospheric water absorption window (1800 nm – 1950 nm) is explicitly sliced out of the dataset to prevent atmospheric scattering from poisoning the covariance matrix calculations.
3. **K-Means Background Clustering:** Instead of treating the landscape as one uniform background, the algorithm segments the scene into distinct spectral clusters (e.g., bare soil, vegetation). 
4. **Adaptive Coherence Estimator (ACE):** For each cluster, a local covariance matrix ($C$) is computed. The ACE algorithm measures the squared Mahalanobis distance between the pixel spectrum ($x$) and the target signature ($d$), normalized by the background variance:
   $$D_{ACE}(x) = \frac{((x - \mu)^T C^{-1} d)^2}{((x - \mu)^T C^{-1} (x - \mu))(d^T C^{-1} d)}$$
   This critical normalization step makes the algorithm **scale-invariant**, effectively immunizing it against false alarms caused by solar illumination changes or bright terrestrial surfaces.

---

## 3. Visual and Quantitative Results

The results generated by the pipeline definitively prove the superiority of the CT-ACE approach. The following visualizations and metrics are auto-generated by the master execution script.

### Stagewise Ablation & Spatial Isolation

The progression of anomaly responses evolves drastically across the evaluated algorithms. 

<p align="center">
  <img src="./updated_output_file/Figure1_Stagewise_Ablation_3.png" width="900" alt="Stagewise Ablation">
</p>
<p align="center">
  <em><b>Figure 1.</b> Stagewise comparison of anomaly responses. Panel (a) shows the original SWIR scene context. Panels (b), (c), and (d) reveal the chaotic noise floors of the CIBR, JRGE, and SFA baselines, where the gas plume is heavily obscured by striping and background artifacts. Panel (e) demonstrates the output of the proposed CT-ACE framework, resulting in a highly cohesive, completely isolated hotspot mask in Panel (f).</em>
</p>

### Profile Analysis & Noise Floor

To quantitatively examine the noise floor of each algorithm, a horizontal 1D spatial cross-section was taken directly through the centroid of the detected CO₂ plume.

<p align="center">
  <img src="./updated_output_file/Figure2_Profile_Analysis_3.png" width="750" alt="Profile Analysis">
</p>
<p align="center">
  <em><b>Figure 2.</b> 1D spatial profile localization. Baseline methods such as JRGE (blue) and SFA (orange) exhibit high turbulence across the entire spatial axis, indicating that they constantly trigger false-positive values over regular terrestrial surfaces. The proposed CT-ACE (Red) maintains a pristine, flat noise floor exactly at 0.0, spiking strictly at the true spatial location of the CO₂ plume.</em>
</p>

### Threshold Sensitivity & False Alarms

A robust algorithm must not rely on a "magic" threshold number to function. The sensitivity of each algorithm was tested across varying percentile thresholds to evaluate false alarm resistance.

<p align="center">
  <img src="./updated_output_file/Figure3_Threshold_Sensitivity_3.png" width="750" alt="Threshold Sensitivity">
</p>
<p align="center">
  <em><b>Figure 3.</b> Algorithm confidence versus background clutter. The baseline CIBR response (grey) is smeared across background noise; meaning that lowering the threshold even slightly triggers a massive volume of false-positive pixels. The proposed CT-ACE algorithm (red) stays nearly flat at zero until the extreme 95th percentile, proving exceptional statistical confidence and resistance to false alarms.</em>
</p>

### Selectivity Metrics (SBR)

While visualizations provide intuitive proof, the **Signal-to-Background Ratio (SBR)** provides rigorous mathematical evidence. SBR calculates how effectively an algorithm amplifies the target gas signature while suppressing background clutter scores.

| Algorithm | Mean Target Score | Mean Background Score | SBR Selectivity |
| :--- | :--- | :--- | :--- |
| CIBR Baseline | 0.10088 | 0.038296 | 2.6343 |
| JRGE Baseline | 0.24659 | 0.137610 | 1.7920 |
| SFA Baseline | 0.10522 | 0.113190 | 0.9296 |
| **CT-ACE (Proposed)** | **0.24416** | **0.005365** | **45.5020** |

*Table 1. Selectivity metrics extracted from the localized plume region versus the total background scene.*

**Conclusion:** Standard baselines achieve an SBR of ~1.0 to 2.6, indicating the detected plume is barely mathematically distinguishable from standard terrain. The proposed CT-ACE framework achieves an SBR of **45.5**, crushing background noise to near absolute zero ($0.005$) and proving its superior background suppression capabilities.

### Geospatial Validation

To fulfill the requirements of practical environmental monitoring, the localized anomalies must be anchored to geographic coordinates to trace emissions back to their terrestrial sources. 

<p align="center">
  <img src="./updated_output_file/Figure4_Geospatial_Validation_3.jpg" width="850" alt="Geospatial Validation">
</p>
<p align="center">
  <em><b>Figure 4.</b> Geographic projection of the thresholded CO₂ hotspot mask. Using the AVIRIS metadata (UTM Zone 11N, WGS-84 datum), the pixel coordinates were mathematically transformed into Latitude and Longitude using the MATLAB Mapping Toolbox. This allows the highly isolated anomaly mask to be overlaid onto a real-world high-resolution satellite basemap for facility identification.</em>
</p>

---

## 4. Toolbox Utilization

This project strictly adheres to MathWorks Challenge guidelines by fully leveraging official, modern MATLAB toolboxes to avoid obsolete, deprecated functions.

* **Hyperspectral Imaging Toolbox:** Data ingestion, wavelength formatting, and spectral band subsetting.
* **Image Processing Toolbox:** Spatial coherence filtering (`medfilt2`) and morphological connected-component analysis (`bwareaopen`).
* **Statistics and Machine Learning Toolbox:** Fast K-means background clustering and Mahalanobis distance covariance matrix calculation.
* **Mapping Toolbox:** Geographic coordinate transformations (`projinv`) and spatial referencing (`maprasterref`, `geoaxes`, `geoscatter`).

---

## 5. Data Acquisition & Setup

### Addressing File Size Constraints
Due to GitHub's strict file size limits (and the impracticality of cloning 41 GB repositories), the raw AVIRIS data file (`f250923t01p00r13_rfl.dat`) is not hosted in this repository. 

To allow reviewers, judges, and researchers to execute the code immediately out-of-the-box, a pre-cropped, memory-safe Region of Interest (ROI) is provided.

### Local Setup Instructions
1. **Clone this repository** to your local machine using `git clone`.
2. Ensure the sample dataset **`proposed_results.mat`** is located in the root directory alongside the `.m` scripts.
3. The codebase utilizes dynamic relative pathing (`pwd`), meaning **no manual directory or path configuration is required.**

---

## 6. Repository Structure & Usage

This project has been heavily refactored to prioritize clean architecture and reproducibility. Fragmented stage files and dummy data have been removed. The entire framework operates via two self-contained master scripts.

### Directory Layout
```text
.
├── generate_all_figures.m       # Master script: Runs entire detection & mapping pipeline
├── histogram_comparison.m       # Statistical validation script: Calculates score distributions
├── proposed_results.mat         # Sample hyperspectral datacube (Required for execution)
├── updated_output_file/         # Auto-generated directory containing all results
│   ├── Figure1_Stagewise_Ablation_3.png
│   ├── Figure2_Profile_Analysis_3.png
│   ├── Figure3_Threshold_Sensitivity_3.png
│   ├── Figure4_Geospatial_Validation_3.jpg
│   └── Table1_Selectivity_Metrics_3.csv
└── README.md                    # Project documentation (this file)
