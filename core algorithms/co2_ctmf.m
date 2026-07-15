function [mf_image, hotspot_mask] = co2_ctmf(mat_file)
% CO2_CTMF  Stage 4 — Cluster-Tuned Matched Filter (CT-ACE, Proposed)
%
% Applies an absorbance-domain transformation and K=4 spectral clustering
% before computing a per-cluster matched filter score (ACE formulation).
% All dependencies are resolved as local functions — this file runs
% standalone without any external helper scripts.
%
% Uses MATLAB Hyperspectral Imaging Toolbox:
%   hypercube()  — data container and wavelength accessor
%
% USAGE
%   [mf_image, hotspot_mask] = co2_ctmf()
%   [mf_image, hotspot_mask] = co2_ctmf(mat_file)
%
% INPUT
%   mat_file  (optional) — full path to proposed_results.mat
%
% OUTPUT
%   mf_image     [rows x cols]  Normalised CT-ACE score map in [0,1]
%   hotspot_mask [rows x cols]  Binary hotspot mask (Otsu, area >= 10 px)
%
% KEY DESIGN CHOICES
%   1. Absorbance transform -log(max(R, 1e-5)) linearises Beer-Lambert
%   2. Water vapour bands 1800–1950 nm excluded from covariance
%   3. K=4 K-Means clusters (rng(1) for reproducibility)
%   4. Per-cluster ACE score: num^2 / (den_bg * den_target)
%   5. Backslash operator (C\d) replaces inv(C)*d throughout
%   6. maprefcells replaces deprecated maprasterref
%
% REFERENCES
%   Marion et al. (2004); Manolakis et al. (2009); Thompson et al. (2015)
% =========================================================================

    clc; close all;

    %% 1. RESOLVE DATA PATH
    if nargin < 1 || isempty(mat_file)
        this_dir = fileparts(mfilename('fullpath'));
        mat_file = fullfile(this_dir, '..', 'proposed_results.mat');
        mat_file = char(java.io.File(mat_file).getCanonicalPath());
    end
    if ~isfile(mat_file)
        error('co2_ctmf: proposed_results.mat not found at:\n  %s', mat_file);
    end

    %% 2. LOAD DATA VIA hypercube()
    fprintf('[CT-ACE] Loading data...\n');
    S = load(mat_file, 'cube', 'wavelengths');
    hcube     = hypercube(double(S.cube), S.wavelengths(:));
    wavelengths = hcube.Wavelength;
    [rows, cols, ~] = size(hcube.DataCube);
    N = rows * cols;

    %% 3. BAND SELECTION — SWIR, exclude water vapour 1800–1950 nm
    valid_bands = (wavelengths >= 1500 & wavelengths <= 2100) & ...
                  ~(wavelengths >= 1800 & wavelengths <= 1950);
    num_valid   = sum(valid_bands);
    wl_valid    = wavelengths(valid_bands);

    %% 4. ABSORBANCE TRANSFORM
    % A(λ) = -log(max(R(λ), 1e-5))
    % Linearises Beer-Lambert: gas signal becomes additive, not multiplicative.
    A_cube   = -log(max(double(hcube.DataCube(:,:,valid_bands)), 1e-5));
    reshapedA = reshape(A_cube, N, num_valid);   % [N x M]

    %% 5. CO2 TARGET SIGNATURE
    co2_sig = buildtargetspectrum(wl_valid);    % [M x 1], L2-normalised

    %% 6. K-MEANS CLUSTERING (K = 4, fixed seed)
    fprintf('[CT-ACE] Clustering (K=4)...\n');
    rng(1);
    [clusterIdx, ~] = kmeans(reshapedA, 4, ...
        'MaxIter',    200, ...
        'Distance',   'sqeuclidean', ...
        'Replicates', 3);

    %% 7. PER-CLUSTER ACE SCORE
    % For cluster k:
    %   Ck = cov(X_k) + 1e-6*I           (regularised covariance)
    %   w  = Ck \ d                       (backslash: avoids explicit inv)
    %   num   = X_c * w                   (matched projection)
    %   den_bg = diag(X_c * (Ck \ X_c')) (per-pixel background energy)
    %   den_t  = d' * w                   (target self-energy)
    %   score  = sign(num) * num^2 / (den_bg * den_t)
    fprintf('[CT-ACE] Computing per-cluster ACE scores...\n');
    ace_output = zeros(N, 1);

    for k = 1:4
        idx  = (clusterIdx == k);
        n_k  = sum(idx);
        if n_k < num_valid + 5
            fprintf('  Cluster %d: only %d pixels — skipped.\n', k, n_k);
            continue;
        end

        X_k  = reshapedA(idx, :);               % [n_k x M]
        mu_k = mean(X_k, 1);
        X_c  = X_k - mu_k;                      % mean-centred

        % Regularised covariance
        Ck   = cov(X_k) + eye(num_valid) * 1e-6;

        % Whitened target vector — BACKSLASH replaces inv(Ck)*co2_sig
        w    = Ck \ co2_sig;                     % [M x 1]

        % ACE numerator: matched filter projection
        num  = X_c * w;                          % [n_k x 1]

        % ACE denominator term 1: per-pixel background energy
        % X_c * Ck^{-1} * X_c^T diagonal
        % Computed as (Ck \ X_c')' to avoid forming full inverse
        X_c_w = (Ck \ X_c')';                   % [n_k x M]
        den1  = sum(X_c .* X_c_w, 2);           % [n_k x 1]

        % ACE denominator term 2: target self-energy
        den2  = co2_sig' * w;                    % scalar

        % ACE score
        ace_output(idx) = sign(num) .* (num.^2) ./ ...
                          (den1 .* den2 + 1e-8);
    end

    %% 8. POST-PROCESSING
    mf_image = reshape(ace_output, rows, cols);
    mf_image = max(mf_image, 0);                 % rectify
    mf_image = medfilt2(mf_image, [3 3]);        % spatial smoothing
    mf_image = mf_image / (max(mf_image(:)) + 1e-10);   % normalise

    %% 9. HOTSPOT MASK — Otsu threshold + area opening
    tau_otsu     = graythresh(mf_image);
    hotspot_mask = bwareaopen(mf_image > tau_otsu, 10);

    % Report metrics
    sig = mean(mf_image(hotspot_mask));
    bg  = mean(mf_image(~hotspot_mask));
    fprintf('[CT-ACE] Coverage: %.2f%%  |  SBR: %.3f  |  Otsu: %.4f\n', ...
        100*sum(hotspot_mask(:))/numel(hotspot_mask), ...
        sig/(bg+1e-10), tau_otsu);

    %% 10. VISUALISATION — score map and hotspot mask
    figure('Name','Stage 4: CT-ACE (Proposed)','Color','w','Position',[160 160 700 340]);
    subplot(1,2,1);
    imagesc(mf_image); axis image; colormap(gca,'parula'); colorbar;
    title('CT-ACE Score Map','FontWeight','bold');

    subplot(1,2,2);
    imshow(hotspot_mask);
    title('Proposed Hotspot Mask (Otsu)','FontWeight','bold');

    %% 11. GEOSPATIAL OVERLAY (called as local function)
    mapgeospatial_overlay(hotspot_mask, rows, cols);
end


% =========================================================================
% LOCAL FUNCTION 1: buildtargetspectrum
% Dual-band Gaussian CO2 template — L2-normalised.
% Kept local so co2_ctmf.m is fully self-contained.
% =========================================================================
function sig = buildtargetspectrum(wavelengths)
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


% =========================================================================
% LOCAL FUNCTION 2: mapgeospatial_overlay
% Projects the binary hotspot mask onto a satellite basemap.
% FIX: replaces deprecated maprasterref with maprefcells.
% =========================================================================
function mapgeospatial_overlay(hotspot_mask, rows, cols)
    % Ground sampling distance and UTM upper-left corner
    dx = 14.4;  dy = 14.4;                     % metres per pixel
    UL_Easting  = 314881.49;
    UL_Northing = 4218923.9;

    % maprefcells — replacement for deprecated maprasterref
    x_limits = [UL_Easting,               UL_Easting  + cols*dx];
    y_limits = [UL_Northing - rows*dy,    UL_Northing           ];
    R = maprefcells(x_limits, y_limits, [rows, cols], ...
                    'ColumnsStartFrom', 'north');

    % Reproject pixel mask to latitude / longitude
    crs       = projcrs(32611);                 % WGS84 UTM Zone 11N
    [X_grid, Y_grid] = worldGrid(R);
    [lat, lon]       = projinv(crs, X_grid, Y_grid);

    % Plot on satellite basemap
    fig = figure('Name','Geospatial Validation','Color','w', ...
                 'Position',[100 100 800 600]);
    gx  = geoaxes;
    geobasemap(gx, 'satellite');
    hold(gx, 'on');
    geoscatter(gx, lat(hotspot_mask), lon(hotspot_mask), ...
               15, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
    geolimits(gx, [min(lat(:)) max(lat(:))], [min(lon(:)) max(lon(:))]);
    title(gx, 'Geospatial Validation — CT-ACE CO_2 Hotspots', ...
          'FontSize',12,'FontWeight','bold');
end
