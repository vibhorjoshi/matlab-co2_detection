function [co2_map, hotspot_mask] = co2_sfa(mat_file)
% CO2_SFA  Stage 3 — Spectral Fitting Algorithm
%
% Measures per-pixel spectral similarity to a dual-band CO2 Gaussian
% absorption template across the SWIR window (1500–2100 nm).
% Uses MATLAB Hyperspectral Imaging Toolbox:
%   hypercube()     — data container and wavelength accessor
%   spectralmatch() — computes Spectral Angle Mapper (SAM) similarity
%
% SBR FIX vs. prior version
%   The original NCC-based score had SBR < 1 because broadband
%   reflectance variation produced positive correlations with the
%   CO2 template across non-plume pixels.  This version uses SAM
%   (inverted and scaled) which measures angular distance in spectral
%   space and is invariant to reflectance magnitude.  A lower SAM angle
%   means closer match to the CO2 template.  The score is converted to
%   a similarity index: score = 1 / (1 + SAM_angle), so high values
%   indicate CO2-like spectra.
%
% USAGE
%   [co2_map, hotspot_mask] = co2_sfa()
%   [co2_map, hotspot_mask] = co2_sfa(mat_file)
%
% INPUT
%   mat_file  (optional) — full path to proposed_results.mat
%
% OUTPUT
%   co2_map      [rows x cols]  Normalised SFA similarity map in [0,1]
%   hotspot_mask [rows x cols]  Binary hotspot mask (Otsu, area >= 10 px)
%
% CO2 TEMPLATE  (buildtargetspectrum local function)
%   Band 1: A=0.30, centre=1575 nm, sigma=15 nm
%   Band 2: A=0.70, centre=2005 nm, sigma=12 nm
%   L2-normalised
%
% REFERENCES
%   Green (2001); Manolakis et al. (2009); Romaniello et al. (2021)
% =========================================================================

    clc; close all;

    %% 1. RESOLVE DATA PATH
    if nargin < 1 || isempty(mat_file)
        this_dir = fileparts(mfilename('fullpath'));
        mat_file = fullfile(this_dir, '..', 'proposed_results.mat');
        mat_file = char(java.io.File(mat_file).getCanonicalPath());
    end
    if ~isfile(mat_file)
        error('co2_sfa: proposed_results.mat not found at:\n  %s', mat_file);
    end

    %% 2. LOAD DATA VIA hypercube()
    fprintf('[SFA] Loading data...\n');
    S = load(mat_file, 'cube', 'wavelengths');
    hcube     = hypercube(double(S.cube), S.wavelengths(:));
    wavelengths = hcube.Wavelength;
    [rows, cols, ~] = size(hcube.DataCube);

    %% 3. BAND SELECTION: SWIR 1500–2100 nm
    swir_mask = wavelengths >= 1500 & wavelengths <= 2100;
    wl_swir   = wavelengths(swir_mask);
    M_swir    = sum(swir_mask);

    %% 4. BUILD CO2 TARGET TEMPLATE
    ref = buildtargetspectrum(wl_swir);     % [M_swir x 1], L2-normalised

    %% 5. SPECTRAL SIMILARITY via spectralmatch() — SAM method
    % spectralmatch returns the Spectral Angle (radians) between each
    % pixel and the reference.  Lower angle = better match.
    fprintf('[SFA] Computing spectralmatch (SAM)...\n');
    swir_cube = double(hcube.DataCube(:,:,swir_mask));   % [rows x cols x M_swir]

    % spectralmatch expects data as [H x W x B] and reference as [1 x B]
    sam_angles = spectralmatch(swir_cube, ref', 'Method', 'sam');
    % sam_angles is [rows x cols], values in [0, pi/2] radians

    % Convert angle to similarity score: 1/(1+angle) in [0,1]
    % Pixels with angle~0 (perfect match) score close to 1.
    % Pixels with angle~pi/2 (orthogonal) score close to 2/pi ≈ 0.64.
    co2_map = 1 ./ (1 + sam_angles);

    %% 6. POST-PROCESSING
    co2_map = medfilt2(co2_map, [3 3]);

    % Normalise to [0,1] using percentile stretch
    P2  = prctile(co2_map(:),  2);
    P98 = prctile(co2_map(:), 98);
    co2_map = (co2_map - P2) / (P98 - P2 + 1e-10);
    co2_map = max(min(co2_map, 1), 0);

    %% 7. HOTSPOT MASK — Otsu threshold + area opening
    tau_otsu     = graythresh(co2_map);
    hotspot_mask = bwareaopen(co2_map > tau_otsu, 10);

    % Report SBR for diagnostic purposes
    sig = mean(co2_map(hotspot_mask));
    bg  = mean(co2_map(~hotspot_mask));
    fprintf('[SFA] Coverage: %.2f%%  |  SBR: %.3f\n', ...
            100*sum(hotspot_mask(:))/numel(hotspot_mask), sig/(bg+1e-10));

    %% 8. VISUALISATION
    figure('Name','Stage 3: SFA','Color','w','Position',[140 140 700 340]);
    subplot(1,2,1);
    imagesc(co2_map); axis image; colormap(gca,'parula'); colorbar;
    title('SFA Spectral Similarity (SAM-based)','FontWeight','bold');

    subplot(1,2,2);
    imshow(hotspot_mask);
    title(sprintf('SFA Hotspots (Otsu = %.3f)',tau_otsu),'FontWeight','bold');
end

% =========================================================================
% LOCAL FUNCTION: buildtargetspectrum
% =========================================================================
function sig = buildtargetspectrum(wavelengths)
% Dual-band Gaussian CO2 absorption template.
% Band 1: amplitude=0.30, centre=1575 nm, sigma=15 nm  (combination band)
% Band 2: amplitude=0.70, centre=2005 nm, sigma=12 nm  (overtone band)
% Output: L2-normalised column vector [length(wavelengths) x 1]
    amp1=0.30; cen1=1575; sig1=15;
    amp2=0.70; cen2=2005; sig2=12;
    sig = zeros(numel(wavelengths), 1);
    for b = 1:numel(wavelengths)
        lam    = wavelengths(b);
        sig(b) = amp1*exp(-0.5*((lam-cen1)/sig1)^2) ...
               + amp2*exp(-0.5*((lam-cen2)/sig2)^2);
    end
    sig = sig / (norm(sig) + 1e-10);
end
