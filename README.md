# 🌍 CO₂ Detection from Hyperspectral Imagery using AVIRIS Data

<p align="center">
  <img src="https://img.shields.io/badge/MATLAB-R2023b%2B-blue.svg" alt="MATLAB">
  <img src="https://img.shields.io/badge/Challenge-MathWorks_Project_251-orange.svg" alt="MathWorks">
  <img src="https://img.shields.io/badge/Data-AVIRIS--Classic-lightblue.svg" alt="AVIRIS">
  <img src="https://img.shields.io/badge/License-MIT-green.svg" alt="License">
  <img src="https://img.shields.io/badge/Status-Refactored_v2-brightgreen.svg" alt="Status">
</p>

> **MathWorks MATLAB & Simulink Challenge — Project #251**  
> *Detection and Visualization of CO₂ Concentration Using Hyperspectral Satellite Data*  
> Authors: Vibhor Joshi · kaveri roy. · dheeraj kumar sharma -IIIT Guwahati  
> 🔗 [MathWorks Challenge Page](https://github.com/mathworks/MATLAB-Simulink-Challenge-Project-Hub/blob/main/projects/Detection%20and%20Visualization%20of%20CO2%20Concentration%20Using%20Hyperspectral%20Satellite%20Data/README.md)

---

## 📑 Table of Contents

1. [Project Overview](#1-project-overview)
2. [Why Standard Filters Fail — The Physical Problem](#2-why-standard-filters-fail--the-physical-problem)
3. [The Progressive Pipeline Architecture](#3-the-progressive-pipeline-architecture)
4. [CT-ACE: The Proposed Algorithm in Detail](#4-ct-ace-the-proposed-algorithm-in-detail)
5. [Visual and Quantitative Results](#5-visual-and-quantitative-results)
6. [MATLAB Toolbox Utilisation](#6-matlab-toolbox-utilisation)
7. [Dataset Acquisition — Step-by-Step](#7-dataset-acquisition--step-by-step)
8. [Repository Structure — Every File Explained](#8-repository-structure--every-file-explained)
9. [Installation and Quick Start](#9-installation-and-quick-start)
10. [Running the Pipeline](#10-running-the-pipeline)
11. [Files to Remove / Not Needed](#11-files-to-remove--not-needed)
12. [Citation and License](#12-citation-and-license)

---

## 1. Project Overview

This repository implements a **Progressive Spectral Conditioning Framework** for detecting atmospheric CO₂ plumes from AVIRIS-Classic Level-2 surface reflectance hyperspectral imagery. The framework addresses a fundamental limitation in standard gas-detection algorithms: atmospheric gases *absorb* specific SWIR wavelengths rather than reflecting them, so standard emissive-target detectors clip the gas signal to zero.

The pipeline cascades four complementary algorithms:

```
AVIRIS Hyperspectral Cube  (224 bands, 365–2496 nm, 14.4 m/px)
          │
          ▼
   ┌─────────────────────────────────────────────────────┐
   │   STAGE 1 — CIBR    Continuum Interpolated Band Ratio│  ← Broad anomaly screening
   └───────────────────────────┬─────────────────────────┘
                               ▼
   ┌─────────────────────────────────────────────────────┐
   │   STAGE 2 — JRGE    Joint Reflectance & Gas Estimator│  ← Background suppression
   └───────────────────────────┬─────────────────────────┘
                               ▼
   ┌─────────────────────────────────────────────────────┐
   │   STAGE 3 — SFA     Spectral Fitting Algorithm       │  ← Dual-band matching
   └───────────────────────────┬─────────────────────────┘
                               ▼
   ┌─────────────────────────────────────────────────────┐
   │   STAGE 4 — CT-ACE  Cluster-Tuned Adaptive Coherence │  ← Final plume extraction
   │                      Estimator                        │
   └───────────────────────────┬─────────────────────────┘
                               ▼
              Final CO₂ Plume Map (Binary Hotspot Mask)
               + Geospatial UTM Overlay + 3D Landscape
```

**Key results on scene `f250923t01p00r13_rfl`:**

| Stage | Algorithm | Coverage | Mean Score | Selectivity |
|:---:|:---|---:|---:|:---:|
| 1 | CIBR | 10.8 % | 0.516 | Moderate |
| 2 | JRGE | 7.5 % | 0.450 | High |
| 3 | SFA | 16.9 % | 0.114 | Low (high recall) |
| 4 | CT-ACE | **0.1 %** | 0.557 | **Very High** |

The cumulative pipeline reduces the anomaly footprint by **99.1 %**, isolating a ~0.2 km² plume-scale region consistent with single-stack industrial emission geometry.

---

## 2. Why Standard Filters Fail — The Physical Problem

### 2.1 The Absorption Sign Problem

Most hyperspectral target detectors were designed for *reflective* targets (minerals, vehicles, rooftops) where pixels of interest have *higher* radiance than the background. Atmospheric CO₂ molecules do the opposite: they *absorb* solar photons at specific SWIR wavelengths, making plume pixels *darker* than their surroundings.

```
Standard matched filter score:
  s = (x - µ)ᵀ Σ⁻¹ d

If the plume makes x darker at absorption bands, the dot product with a
positive template d returns a NEGATIVE score. Researchers clip s < 0 to
zero, which erases the gas signature entirely.

CT-ACE solution: Convert to absorbance space first:
  A = −log(R)
  
In absorbance, absorption events become positive additive spikes.
```

### 2.2 The Scale-Variance Problem

A global matched filter estimates a single covariance matrix Σ from the entire scene. Over heterogeneous terrain (bare rock, vegetation, shadow, concrete), the single Σ overestimates background variance. A plume pixel in a homogeneous sub-region produces a small detection score relative to the inflated Σ, causing missed detections and false alarms.

### 2.3 The JRGE Zero-Output Bug (Fixed in This Repository)

The original JRGE implementation passed wavelengths in **nanometres** to MATLAB's `csaps()` function. The smoothing parameter `p` is scale-sensitive:

```
Transition region:  p ≈ 1 / (1 + ε)   where ε = h³/16

For h ≈ 10 nm  →  ε ≈ 62.5  →  transition at p ≈ 0.016
Setting p = 0.90 in nm-space → spline acts as near-interpolant (p→1)
The spline fits THROUGH the CO₂ absorption dip → residual = ZERO
```

**Fix applied in `src/co2_jrge_fixed.m`:** wavelengths converted to **micrometres** before `csaps()`. At h ≈ 0.010 µm, ε ≈ 6.25×10⁻⁸, the transition shifts to p ≈ 0.9999, and p = 0.90 correctly produces a smooth spline that spans the absorption feature without fitting into it.

---

## 3. The Progressive Pipeline Architecture

Each stage reduces a different class of false detection. The cascade is designed so that each algorithm operates on **increasingly conditioned** spectral data.

### Stage 1 — CIBR: Continuum Interpolated Band Ratio

**Purpose:** Rapid, low-cost initial screening of pixels exhibiting absorption-like behaviour near 2.05 µm.

**Physical basis:** In the absence of elevated CO₂, surface reflectance varies smoothly across the 2.0–2.1 µm range. A gas plume causes a localised dip below the linearly interpolated spectral continuum.

**Equations:**
```
Continuum:    C(r,c) = [ R(r,c, λ_L) + R(r,c, λ_R) ] / 2

              λ_L = 2010 nm  (left shoulder)
              λ_A = 2050 nm  (CO₂ absorption centre)
              λ_R = 2090 nm  (right shoulder)

CIBR score:   CIBR_raw(r,c) = C(r,c) − R(r,c, λ_A)

Normalised:   CIBR(r,c) = (CIBR_raw − P₂) / (P₉₈ − P₂ + ε)
```

**Why difference, not ratio?** At low-reflectance pixels (shadow, water), `R_A ≈ 0` makes the ratio `C/R_A` numerically unstable. The difference formulation is bounded by `C` and never diverges.

**Output:** High-recall binary hotspot mask (~10.8 % coverage). Intentionally broad — designed for maximum recall, not precision.

**Limitation addressed at next stage:** Carbonate minerals and atmospheric water vapour produce absorption near 2.0 µm, generating false positives that pass the Otsu threshold.

---

### Stage 2 — JRGE: Joint Reflectance and Gas Estimator

**Purpose:** Decompose the observed spectrum into a smooth surface reflectance component and a narrowband gas residual, suppressing mineralogical false positives.

**Physical basis:** Surface reflectance varies gradually over tens of nanometres. CO₂ absorption is a spectrally narrow feature (≈20–50 nm). A smooth spline fitted to the full-band spectrum captures the broadband surface component; the residual in the CO₂ window isolates the gas signal.

**Equations:**
```
Initial spline:   r̂ᵢ⁽⁰⁾ = S₀.₉₀[rᵢ, λ_µm]     (λ MUST be in µm)

Gas estimate:     ĝᵢ⁽⁰⁾ = median(r̂ᵢ|_ℬ) − median(rᵢ|_ℬ)
                  ℬ = {b : λ_b ∈ [2000, 2100] nm}

Iterative refine: r̂ᵢ⁽ᵗ⁺¹⁾ = r̂ᵢ⁽ᵗ⁾ + α(rᵢ − r̂ᵢ⁽ᵗ⁾),  α=0.05, T=3

Final map:        G(r,c) = (ĝ⁽³⁾ − g_min) / (g_max + ε)
```

**Key parameter choices:**
- `p = 0.90` (in µm): smooth enough to span the narrow CO₂ dip without fitting into it
- `α = 0.05`: conservative step prevents over-correction that would absorb the gas signal
- `T = 3`: cumulative update ≤14.3 % of initial residual; prevents convergence to zero

**Output:** Gas density map + background-subtracted residual cube passed to CT-ACE as the conditioned input.

---

### Stage 3 — SFA: Spectral Fitting Algorithm

**Purpose:** Extend detection from three bands (CIBR) to the full 1500–2100 nm SWIR window, exploiting both CO₂ absorption features simultaneously.

**Physical basis:** CO₂ absorbs at two distinct SWIR bands:
- **1.575 µm** — 3ν₁+ν₃ combination mode (weaker, A₁=0.30)
- **2.005 µm** — 4ν₃+2δ combination mode (dominant, A₂=0.70)

Correlated depressions at *both* wavelengths simultaneously is a stronger indicator of CO₂ than a single-band dip (lower probability of mineralogical coincidence).

**Equations:**
```
Dual-Gaussian template:
d(λ) = 0.30·exp[−(λ−1575)²/(2·15²)] + 0.70·exp[−(λ−2005)²/(2·12²)]
d ← d / ‖d‖   (unit normalise)

Pixel standardisation: x̃ = (x − µ_x) / σ_x

NCC score: s(r,c) = x̃ · (d − d̄) / N_SWIR
```

**Relationship to matched filter:** SFA is mathematically equivalent to the linearised MF with Σ = I (identity covariance). The identity assumption ignores spectral correlations, reducing computation from O(B³) to O(B), at the cost of not whitening structured background noise.

**Why coverage increases to 16.9 %:** Water vapour near 1.9 µm partially correlates with the d(λ₁) component. This is intentional — SFA is designed as a **high-recall** stage that provides a broad candidate pool for CT-ACE.

---

### Stage 4 — CT-ACE: Cluster-Tuned Adaptive Coherence Estimator

See [Section 4](#4-ct-ace-the-proposed-algorithm-in-detail) for the full derivation.

---

## 4. CT-ACE: The Proposed Algorithm in Detail

### 4.1 Why CTMF is Not Enough

The original Cluster-Tuned Matched Filter (Marion et al., 2004) partitions the scene and applies a local matched filter per cluster:

```
s_i = (x_i − µ_k)ᵀ Cₖ⁻¹ d
```

This is an unnormalised projection. The score magnitude depends on:
- The energy of the pixel deviation `(x_i − µ_k)`
- The energy of the target vector `d`
- Both of which vary across clusters and illumination conditions

Bright surfaces (concrete, salt flats) produce large `(x − µ)` vectors and therefore large CTMF scores *regardless* of spectral similarity to CO₂, causing false alarms.

### 4.2 The ACE Score: Scale-Invariant Detection

The Adaptive Coherence Estimator normalises the CTMF score by both the pixel background energy and the target energy:

```
                   [ (x − µ_k)ᵀ Cₖ⁻¹ d ]²
D_ACE(x) = ─────────────────────────────────────────────────
            [(x − µ_k)ᵀ Cₖ⁻¹ (x − µ_k)] · [dᵀ Cₖ⁻¹ d]
```

**Physical interpretation:**
- Numerator: squared projection of background-whitened pixel onto target direction
- First denominator term: total background-whitened pixel energy
- Second denominator term: target self-energy under cluster covariance
- Result: bounded [0, 1], where 1 = perfect spectral match to CO₂ template

**Why this eliminates false alarms:** A bright concrete pixel has high `(x−µ)ᵀCₖ⁻¹(x−µ)` (large background energy) as well as high numerator energy. The ratio remains small unless the spectral *shape* of the deviation matches the CO₂ template. The algorithm is immunised against solar illumination gradients and surface albedo variation.

### 4.3 Beer-Lambert Pre-Transform

Before clustering and filtering, reflectance is converted to absorbance:

```
A(r, c, λ) = −log₁₀[ R(r, c, λ) + ε ],   ε = 10⁻⁶
```

In absorbance space, CO₂ absorption at 2.05 µm appears as a **positive additive spike** rather than a negative dip. The matched filter template `d` therefore has the correct sign, and no negative clipping is needed.

### 4.4 Water Vapour Band Exclusion

The opaque atmospheric water vapour window (1800–1950 nm) contains zero useful signal but contributes high-variance noise to the covariance matrix estimation. These bands are excluded before k-means clustering:

```matlab
validBands = wavelength < 1800 | wavelength > 1950;
cube_filtered = absorbanceCube(:,:,validBands);
```

### 4.5 Cluster Covariance with Tikhonov Regularisation

```
Cluster partitioning:  c_i = argmin_k ‖x_i − µ_k‖²,   K=4, seed=1

Regularised covariance: Cₖ = Cov(X_k) + λI,   λ = 10⁻⁶

MF weight vector:       w_k = Cₖ⁻¹ d   (solved via Cₖ\d, no explicit inverse)

CT-ACE score:           D_ACE computed per cluster, assembled into full map

Final mask:             M = S > τ,   τ chosen by Otsu or P95
```

**Why K=4?** AVIRIS scenes over mixed terrain typically contain 4–5 spectrally dominant surface types (photosynthetically active vegetation, dry soil/rock, dark surfaces, bright impervious cover). K too small conflates distinct surfaces; K too large creates under-populated clusters where Cov(Xₖ) is rank-deficient.

**Why λ=10⁻⁶?** At AVIRIS-Classic 10 nm sampling and 200×200 working array, each cluster has ~10,000 pixels against 224 bands. The sample covariance is well-conditioned; λ=10⁻⁶ is several orders of magnitude below typical diagonal entries (~10⁻³) and provides numerical stability without distorting the filter.

### 4.6 CT-ACE vs Standard CTMF vs Global MF — Summary

| Property | Global MF | CTMF (Marion 2004) | **CT-ACE (Proposed)** |
|:---|:---:|:---:|:---:|
| Background model | Global Σ | Per-cluster Cₖ | Per-cluster Cₖ |
| Score normalisation | None | None | By pixel + target energy |
| Scale invariance | ✗ | ✗ | **✓** |
| Absorbance transform | ✗ | ✗ | **✓** |
| Water vapour exclusion | ✗ | ✗ | **✓** |
| False alarm resistance | Low | Moderate | **High** |
| MODTRAN required | Optional | Optional | **Not required** |
| Final coverage | ~10 % | ~5 % | **~0.1 %** |


## 6. MATLAB Toolbox Utilisation

| Toolbox | Functions Used | Purpose |
|:---|:---|:---|
| **Hyperspectral Imaging Library** | `hypercube()` | Standardised AVIRIS data ingestion (replaces custom `multibandread` reader) |
| **Image Processing Toolbox** | `medfilt2`, `bwareaopen`, `graythresh`, `bwconncomp`, `regionprops`, `label2rgb` | Spatial filtering, morphological analysis, Otsu thresholding |
| **Statistics & ML Toolbox** | `kmeans`, `cov`, `prctile` | K-means background clustering, covariance estimation, percentile thresholds |
| **Curve Fitting Toolbox** | `csaps`, `fnval` | Cubic smoothing spline for JRGE background estimation |
| **Mapping Toolbox** | `projinv`, `projcrs`, `maprasterref`, `mapshow`, `geoshow` | UTM→lat/lon reprojection (`projinv` replaces removed `minvtran`), geospatial overlay |

---

## 7. Dataset Acquisition — Step-by-Step

### About the Dataset

**Scene:** `f250923t01p00r13_rfl`  
**Format:** ENVI Band Interleaved by Line (BIL), 32-bit IEEE float  
**Size:** 1937 samples × 24068 lines × 224 bands (~41.8 GB uncompressed)  
**Resolution:** 14.4 m/pixel, UTM Zone 11N, WGS-84  
**UL Corner:** Easting 577,561.59 m, Northing 4,228,899.2 m  
**Wavelengths:** 378.93 – 2498.34 nm (224 bands, ~10 nm sampling)  

### Option A — Quick Start with Pre-computed Results (Recommended)

If you do not have access to the 41.8 GB AVIRIS file, a pre-computed `.mat` file is provided containing all four stage outputs for a 200×200 working array.

```bash
# 1. Clone the repository
git clone https://github.com/vibhorjoshi/co2-detection-matlab.git
cd co2-detection-matlab

# 2. The proposed_results.mat is already in /data
#    Run validation directly
```

```matlab
% In MATLAB
cd co2-detection-matlab
entry_point_verification   % auto-detects proposed_results.mat and runs
```

### Option B — Full Pipeline with Raw AVIRIS Data

**Step 1: Download from JPL AVIRIS Data Portal**
1. Navigate to https://aviris.jpl.nasa.gov/dataportal/
2. Search for flight line `f250923t01p00r13`
3. Select the **Level-2 Surface Reflectance** product
4. Download both files:
   - `f250923t01p00r13_rfl` (binary, ~41.8 GB)
   - `f250923t01p00r13_rfl.hdr` (plain text header, ~8 KB)

**Step 2: Organise files**
```
co2-detection-matlab/
└── data/
    ├── f250923t01p00r13_rfl          ← binary image (41.8 GB)
    ├── f250923t01p00r13_rfl.hdr      ← ENVI header
    └── proposed_results.mat          ← pre-computed (already provided)
```

**Step 3: Verify header metadata**

Open `f250923t01p00r13_rfl.hdr` in any text editor and confirm:
```
samples   = 1937
lines     = 24068
bands     = 224
interleave = bil
data type = 4      (32-bit float)
map info  = {UTM, 1, 1, 577561.590, 4228899.200, 14.400, 14.400, 11, North, WGS-84}
```

**Step 4: Run**
```matlab
cd co2-detection-matlab
entry_point_verification    % detects binary file → runs full pipeline
```

> ⚠️ **Memory note:** The full 41.8 GB scene is never fully loaded. The pipeline crops 1000×1000 pixels via `multibandread` range selectors, requiring approximately 900 MB RAM for the raw crop, which is reduced to ~85 MB after downsampling ×5.

### Option C — USGS EarthExplorer (Alternative)

1. Visit https://earthexplorer.usgs.gov
2. Create a free account
3. Under **Data Sets → Aerial Imagery → AVIRIS**, search for the scene using the flight line ID

---

## 8. Repository Structure — Every File Explained

```
co2-detection-matlab/
│
├── 📄 entry_point_verification.m          ← MASTER SCRIPT — run this first
│
├── 📁 config/
│   └── pipeline_config.m                  ← All parameters, paths, constants
│
├── 📁 src/                                ← Core algorithm implementations
│   ├── buildtargetspectrum.m              ← Dual-Gaussian CO₂ template generator
│   ├── co2_jrge_fixed.m                   ← JRGE (critical bug fix: nm→µm)
│   ├── computectmf.m                      ← CT-ACE / CTMF (fixed fn name)
│   ├── mapgeospatial_overlay.m            ← Geospatial projection (projinv)
│   ├── profileanalysis.m                  ← 1-D horizontal score profile
│   └── sixstage_pipeline.m               ← Master cascade orchestrator
│
├── 📁 validation/                         ← Comparison and validation scripts
│   ├── baseline_ctmf.m                    ← Baseline: raw cube → CTMF only
│   ├── save_proposed_results.m            ← Run all 4 stages, save .mat
│   ├── figure_comparison.m               ← 6-panel comparison figure
│   ├── histogram_comparison.m            ← Score distribution histograms
│   ├── compute_selectivity_metrics.m     ← CSV quantitative metrics table
│   ├── ablation_study.m                  ← 4-case ablation with CSV output
│   ├── ablation_figure.m                 ← 2×2 ablation map panels
│   ├── threshold_sensitivity.m           ← Otsu/P85/P90/P95 comparison
│   └── geospatial_validation.m           ← UTM-referenced figures + 3D
│
├── 📁 legacy/                             ← Original submissions (kept for reference)
│   ├── co2_cibr.m                         ← Original CIBR (hardcoded paths)
│   ├── co2_jrge.m                         ← Original JRGE (zero-output bug)
│   ├── co2_sfa.m                          ← Original SFA (has clc/clear bug)
│   └── main_co2_visualisation.m          ← Original visualisation script
│
├── 📁 data/
│   ├── proposed_results.mat              ← Pre-computed results (10 MB, provided)
│   ├── f250923t01p00r13_rfl              ← Raw AVIRIS binary (not on GitHub)
│   └── f250923t01p00r13_rfl.hdr         ← ENVI header (included)
│
├── 📁 results/                            ← Auto-generated outputs
│   ├── geo_fig01_rgb_overview.png
│   ├── geo_fig02_swir_false_colour.png
│   ├── geo_fig03_pipeline_progression.png
│   ├── geo_fig04_ctmf_heatmap.png
│   ├── geo_fig05_connected_components.png
│   ├── geo_fig06_3d_landscape.png
│   ├── geo_fig07_threshold_sensitivity.png
│   ├── geo_fig08_composite_validation.png
│   ├── comparison_pipeline.png
│   ├── histogram_scores.png
│   ├── ablation_maps.png
│   ├── quantitative_metrics.csv
│   ├── ablation_results.csv
│   └── threshold_results.csv
│
├── 📁 docs/
│   └── co2_detection.pdf                 ← Full manuscript (IEEE format)
│
└── 📄 LICENSE                            ← MIT License
```

### File-by-File Purpose

| File | Role | Key Detail |
|:---|:---|:---|
| `entry_point_verification.m` | Single entry point | Auto-detects run mode (full/precomputed/demo); no path editing needed |
| `config/pipeline_config.m` | Central config | All 30+ parameters in one place; paths relative to repo root |
| `src/buildtargetspectrum.m` | CO₂ template | Dual-Gaussian with HITRAN parameters; validates SWIR window coverage |
| `src/co2_jrge_fixed.m` | JRGE (fixed) | Converts λ to µm before `csaps`; outputs residual cube for CT-ACE |
| `src/computectmf.m` | CT-ACE / CTMF | Function name matches filename (fixed mismatch); targets CO₂ at 2.0 µm |
| `src/mapgeospatial_overlay.m` | Geo overlay | Uses `projinv` (R2022b+); fallback to UTM-axis display if no Mapping Toolbox |
| `src/profileanalysis.m` | 1-D profile | Auto-selects peak row; computes PBR, half-peak width, contrast |
| `src/sixstage_pipeline.m` | Master pipeline | True progressive cascade (JRGE residuals fed to CT-ACE, not raw cube) |
| `validation/baseline_ctmf.m` | Baseline | Applies CTMF to raw reflectance with single-band Gaussian target |
| `validation/save_proposed_results.m` | Proposed run | Saves all 4 intermediate score maps to `.mat` |
| `validation/geospatial_validation.m` | Geo validation | Loads BIL file, runs pipeline, generates 8 georeferenced figures |
| `data/f250923t01p00r13_rfl.hdr` | AVIRIS metadata | Included; needed even when binary not available (for wavelength vector) |
| `data/proposed_results.mat` | Pre-computed | Enables validation scripts without the 41.8 GB binary |

---

## 9. Installation and Quick Start

### Requirements

| Requirement | Version |
|:---|:---|
| MATLAB | R2023b or later |
| Image Processing Toolbox | Any recent |
| Statistics and Machine Learning Toolbox | Any recent |
| Curve Fitting Toolbox | Any recent |
| Hyperspectral Imaging Library | R2020b+ |
| Mapping Toolbox | R2022b+ (for `projinv`) |

> Mapping Toolbox is **optional**: geospatial figures fall back to UTM-axis display if not licensed.

### Installation

```bash
# Clone the repository
git clone https://github.com/vibhorjoshi/co2-detection-matlab.git
cd co2-detection-matlab
```

No toolbox installation beyond standard MATLAB is required. All `addpath` calls are handled automatically by `entry_point_verification.m`.

### Verify installation

```matlab
% From the repository root directory
entry_point_verification
% Expected output (demo mode if no data file):
%   Repository root : /path/to/co2-detection-matlab
%   Execution mode  : DEMO
%   Synthetic demo complete.
```

---

## 10. Running the Pipeline

### Full run (with AVIRIS binary)
```matlab
entry_point_verification
% Automatically:
%   1. Adds src/ and config/ to path
%   2. Detects binary + HDR → runs sixstage_pipeline
%   3. Saves all stage outputs to results/proposed_results.mat
%   4. Generates all 8 geospatial figures
%   5. Runs profileanalysis on peak row
%   6. Runs mapgeospatial_overlay
```

### Run only validation (pre-computed data)
```matlab
% Requires: data/proposed_results.mat
addpath('src'); addpath('config');

% Score distributions
cd validation
histogram_comparison

% Ablation study
ablation_study
ablation_figure

% Threshold sensitivity
threshold_sensitivity

% Geospatial figures (requires HDR file)
geospatial_validation
```

### Run individual stages
```matlab
addpath('src'); addpath('config');
cfg = pipeline_config();

% Load a crop manually
hc  = hypercube(fullfile(cfg.Path.Data, cfg.Path.DataFile), ...
                fullfile(cfg.Path.Data, cfg.Path.HeaderFile));
cube = double(hc.DataCube(1:200, 1:200, :)) / cfg.Sensor.ScaleFactor;
wl   = double(hc.Wavelength);

% Stage 2: JRGE (fixed)
cfg.wavelength = wl;
[jrgeScore, jrgeMask, residualCube] = co2_jrge_fixed(cube, cfg);

% Stage 4: CT-ACE on residuals
[d, ~, ~] = buildtargetspectrum(wl, cfg);
[S, M]    = computectmf(residualCube, d(1:size(residualCube,3)), ...
                        cfg.Params.CTMF.Clusters, cfg.Params.CTMF.Regularization);
```

---

---

## 11. Citation and License

### Citing This Work

```bibtex
@misc{joshi2025co2aviris,
  title   = {A Multi-Stage Hyperspectral Framework for CO$_2$ Plume Detection 
             Using CIBR, JRGE, SFA, and CT-ACE from AVIRIS Imagery},
  author  = {Joshi, Vibhor and Y., Namratha},
  year    = {2025},
  note    = {MathWorks MATLAB Simulink Challenge Project \#251},
  url     = {https://github.com/vibhorjoshi/co2-detection-matlab}
}
```

### Key References

| Reference | Role in this work |
|:---|:---|
| Marion et al. (2004), *IEEE TGRS* | Original CTMF algorithm baseline |
| Kim et al. (2025), *JGR Atmospheres* | Lognormal MF benchmark; POD analysis |
| Romaniello et al. (2021), *Remote Sens.* | CIBR methodology and PRISMA validation |
| Dennison et al. (2013), *Remote Sens. Environ.* | First facility-scale CO₂ AVIRIS mapping |
| Otsu (1979), *IEEE Trans. SMC* | Thresholding methodology |
| HITRAN 2020 Database | CO₂ absorption parameters (λ₁, λ₂, σ₁, σ₂) |

### License

This project is released under the **MIT License**.  
Copyright © 2025 Vibhor Joshi, Namratha Y., IIIT Guwahati.  
See `LICENSE` for full terms.

---

<p align="center">
  <sub>Built for MathWorks MATLAB & Simulink Challenge Project #251 · IIIT Guwahati · 2025</sub>
</p>
