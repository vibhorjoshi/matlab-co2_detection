function [mf_image, hotspot_mask] = computectmf(dataFile, headerFile)
% COMPUTECTMF  Stage 4 — Cluster-Tuned Matched Filter (CT-ACE, Proposed)
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
%   [mf_image, hotspot_mask] = computectmf(dataFile, headerFile)
%
% INPUT
%   dataFile     — full path to the raw AVIRIS binary file (.dat / .bin)
%   headerFile   — full path to the AVIRIS header file (.hdr)
%
% OUTPUT
%   mf_image     [rows x cols]  Normalised CT-ACE score map in [0,1]
%   hotspot_mask [rows x cols]  Binary hotspot mask (Otsu, area >= 10 px)
% =========================================================================

    clc; close all;

    %% 1. LOAD DATA VIA hypercube()
    if nargin < 2
        error('Please provide both the data file and header file paths.');
    end
    
    fprintf('[CT-ACE] Loading raw data using Hyperspectral Imaging Library...\n');
    hcube = hypercube(dataFile, headerFile);
    wavelengths = hcube.Wavelength;
    
    % Use a cropped region for memory-safe processing (optional, adjust as needed)
    cube_data = double(hcube.DataCube(1:200, 1:200, :));
    [rows, cols, ~] = size(cube_data);
    N = rows * cols;

    %% 2. BAND SELECTION — SWIR, exclude water vapour 1800–1950 nm
    valid_bands = (wavelengths >= 1500 & wavelengths <= 2100) & ...
                  ~(wavelengths >= 1800 & wavelengths <= 1950);
    num_valid   = sum(valid_bands);
    wl_valid    = wavelengths(valid_bands);

    %% 3. ABSORBANCE TRANSFORM
    % A(λ) = -log(max(R(λ), 1e-5))
    % Linearises Beer-Lambert: gas signal becomes additive, not multiplicative.
    A_cube = -log(max(cube_data(:,:,valid_bands), 1e-5));
    reshapedA = reshape(A_cube, N, num_valid);   % [N x M]

    %% 4. CO2 TARGET SIGNATURE
    co2_sig = buildtargetspectrum(wl_valid);    % [M x 1], L2-normalised

    %% 5. K-MEANS CLUSTERING (K = 4, fixed seed)
    fprintf('[CT-ACE] Clustering (K=4)...\n');
    rng(1);
    [clusterIdx, ~] = kmeans(reshapedA, 4, ...
        'MaxIter',    200, ...
        'Distance',   'sqeuclidean', ...
        'Replicates', 3);

    %% 6. PER-CLUSTER ACE SCORE
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
        % Computed as (Ck \ X_c')' to avoid forming full inverse matrix
        X_c_w = (Ck \ X_c')';                   % [n_k x M]
        den1  = sum(X_c .* X_c_w, 2);           % [n_k x 1]

        % ACE denominator term 2: target self-energy
        den2  = co2_sig' * w;                    % scalar

        % ACE score
        ace_output(idx) = sign(num) .* (num.^2) ./ (den1 .* den2 + 1e-8);
    end

    %% 7. POST-PROCESSING
    mf_image = reshape(ace_output, rows, cols);
    mf_image = max(mf_image, 0);                 % rectify
    mf_image = medfilt2(mf_image, [3 3]);        % spatial smoothing
    mf_image = mf_image / (max(mf_image(:)) + 1e-10);   % normalise

    %% 8. HOTSPOT MASK — Otsu threshold + area opening
    tau_otsu     = graythresh(mf_image);
    hotspot_mask = bwareaopen(mf_image > tau_otsu, 10);

    % Report metrics
    sig = mean(mf_image(hotspot_mask));
    bg  = mean(mf_image(~hotspot_mask));
    fprintf('[CT-ACE] Coverage: %.2f%%  |  SBR: %.3f  |  Otsu: %.4f\n', ...
        100*sum(hotspot_mask(:))/numel(hotspot_mask), ...
        sig/(bg+1e-10), tau_otsu);

    %% 9. VISUALISATION
    figure('Name','Stage 4: CT-ACE (Proposed)','Color','w','Position',[160 160 700 340]);
    subplot(1,2,1);
    imagesc(mf_image); axis image; colormap(gca,'parula'); colorbar;
    title('CT-ACE Score Map','FontWeight','bold');

    subplot(1,2,2);
    imshow(hotspot_mask);
    title('Proposed Hotspot Mask (Otsu)','FontWeight','bold');

    %% 10. GEOSPATIAL OVERLAY (called as local function)
    mapgeospatial_overlay(hotspot_mask, rows, cols);
end


% =========================================================================
% LOCAL FUNCTION 1: buildtargetspectrum
% Dual-band Gaussian CO2 template — L2-normalised.
% Kept local so computectmf.m is fully self-contained.
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
% FIX: Space removed from MarkerFaceAlpha property.
% =========================================================================
function mapgeospatial_overlay(hotspot_mask, rows, cols)
    % Ground sampling distance and UTM upper-left corner
    dx = 14.4;  dy = 14.4;                     % metres per pixel
    UL_Easting  = 577561.590;                  % AVIRIS f250923t01p00r13 Easting
    UL_Northing = 4228899.200;                 % AVIRIS f250923t01p00r13 Northing

    % maprefcells — replacement for deprecated maprasterref
    x_limits = [UL_Easting,               UL_Easting  + cols*dx];
    y_limits = [UL_Northing - rows*dy,    UL_Northing           ];
    R = maprefcells(x_limits, y_limits, [rows, cols], ...
                    'ColumnsStartFrom', 'north');

    % Reproject pixel mask to latitude / longitude
    crs = projcrs(32611);                 % WGS84 UTM Zone 11N
    [X_grid, Y_grid] = worldGrid(R);
    [lat, lon] = projinv(crs, X_grid, Y_grid);

    % Plot on satellite basemap
    fig = figure('Name','Geospatial Validation','Color','w', ...
                 'Position',[100 100 800 600]);
    gx  = geoaxes;
    geobasemap(gx, 'satellite');
    hold(gx, 'on');
    
    % Syntactical Fix Applied: MarkerFaceAlpha is a single unspaced keyword
    geoscatter(gx, lat(hotspot_mask), lon(hotspot_mask), ...
               15, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
               
    geolimits(gx, [min(lat(:)) max(lat(:))], [min(lon(:)) max(lon(:))]);
    title(gx, 'Geospatial Validation — CT-ACE CO_2 Hotspots', ...
          'FontSize',12,'FontWeight','bold');
end
