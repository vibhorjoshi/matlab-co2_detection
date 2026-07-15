function [co2_index, hotspot_mask] = co2_cibr(mat_file)
% CO2_CIBR  Stage 1 — Continuum Interpolated Band Ratio
%
% Detects CO2 absorption at ~2.05 µm by measuring the depth of the
% observed reflectance below a locally linear spectral continuum.
% Uses MATLAB Hyperspectral Imaging Toolbox:
%   hypercube()       — data container and wavelength accessor
%   removecontinuum() — convex-hull continuum normalisation
%
% USAGE
%   [co2_index, hotspot_mask] = co2_cibr()
%   [co2_index, hotspot_mask] = co2_cibr(mat_file)
%
% INPUT
%   mat_file  (optional) — full path to proposed_results.mat
%             Default: <repo_root>/proposed_results.mat
%
% OUTPUT
%   co2_index    [rows x cols]  Normalised CIBR score map in [0,1]
%   hotspot_mask [rows x cols]  Binary hotspot mask (Otsu threshold,
%                               components < 10 px removed)
%
% BAND CENTRES (AVIRIS-Classic alignment, avoids 1.9 µm water band)
%   Left shoulder  : 2010 nm
%   Absorption peak: 2050 nm
%   Right shoulder : 2090 nm
%
% REFERENCES
%   Green (2001); Romaniello et al. (2021)
% =========================================================================

    clc; close all;

    %% 1. RESOLVE DATA PATH
    if nargin < 1 || isempty(mat_file)
        % Locate repo root two levels above this file
        this_dir = fileparts(mfilename('fullpath'));
        mat_file = fullfile(this_dir, '..', 'proposed_results.mat');
        mat_file = char(java.io.File(mat_file).getCanonicalPath());
    end
    if ~isfile(mat_file)
        error('co2_cibr: proposed_results.mat not found at:\n  %s', mat_file);
    end

    %% 2. LOAD DATA VIA hypercube()
    fprintf('[CIBR] Loading data...\n');
    S = load(mat_file, 'cube', 'wavelengths');
    % Wrap in hypercube object — enables Toolbox functions
    hcube = hypercube(double(S.cube), S.wavelengths(:));
    wavelengths = hcube.Wavelength;
    [rows, cols, ~] = size(hcube.DataCube);

    %% 3. CONTINUUM REMOVAL via removecontinuum()
    % removecontinuum normalises each pixel spectrum by its convex-hull
    % continuum.  CR(λ) = R(λ) / Continuum(λ).
    % At an absorption band, CR < 1; depth = 1 - CR gives the band depth.
    fprintf('[CIBR] Computing continuum removal...\n');
    hcube_cr = removecontinuum(hcube);          % [rows x cols x bands], values in [0,1]
    cr_cube   = hcube_cr.DataCube;

    %% 4. BAND SELECTION (nearest available wavelength)
    [~, b_A] = min(abs(wavelengths - 2050));    % absorption peak
    % CIBR score = continuum-removal band depth at the absorption centre
    % (1 - CR(lambda_A)) is positive where CO2 absorbs and zero elsewhere
    cibr_raw = 1 - cr_cube(:,:,b_A);

    %% 5. MANUAL THREE-BAND SCORE (fallback / cross-check)
    % Traditional CIBR: CIBR = ((R_L + R_R)/2) - R_A
    % Kept as a cross-check against the continuum-removal approach above.
    [~, b_L] = min(abs(wavelengths - 2010));
    [~, b_R] = min(abs(wavelengths - 2090));
    R = double(hcube.DataCube);
    cibr_3band = ((max(R(:,:,b_L),0) + max(R(:,:,b_R),0)) / 2) ...
                 - max(R(:,:,b_A), 0);
    cibr_3band = max(cibr_3band, 0);

    % Use continuum-removal score as the primary output;
    % ensure non-negative
    co2_index = max(cibr_raw, 0);

    %% 6. POST-PROCESSING
    co2_index = medfilt2(co2_index, [3 3]);

    % Percentile-based normalisation (robust to outliers)
    P2  = prctile(co2_index(:),  2);
    P98 = prctile(co2_index(:), 98);
    co2_index = (co2_index - P2) / (P98 - P2 + 1e-10);
    co2_index = max(min(co2_index, 1), 0);

    %% 7. HOTSPOT MASK — Otsu threshold + area opening
    tau          = graythresh(co2_index);           % Otsu
    hotspot_mask = bwareaopen(co2_index > tau, 10);

    fprintf('[CIBR] Coverage: %.2f%%  |  Otsu threshold: %.4f\n', ...
            100*sum(hotspot_mask(:))/numel(hotspot_mask), tau);

    %% 8. VISUALISATION
    figure('Name','Stage 1: CIBR','Color','w','Position',[100 100 900 380]);
    subplot(1,3,1);
    imagesc(R(:,:,b_A)); axis image; colormap(gca,'gray');
    title(sprintf('Scene at %.0f nm',wavelengths(b_A)),'FontWeight','bold');

    subplot(1,3,2);
    imagesc(co2_index); axis image; colormap(gca,'parula'); colorbar;
    title('CIBR Score (normalised)','FontWeight','bold');

    subplot(1,3,3);
    imshow(hotspot_mask);
    title(sprintf('CIBR Hotspots (Otsu = %.3f)',tau),'FontWeight','bold');
end
