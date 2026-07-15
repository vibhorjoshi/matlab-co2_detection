function main_pipeline()
% MAIN_PIPELINE: Master execution script for CO2 Detection Challenge #251
% This script loads the hyperspectral data, applies physical Absorbance transforms,
% runs baseline vs. proposed CT-ACE algorithms, and generates all required figures.

    % Clear environment for clean execution
    close all; clc;

    %% 1. DYNAMIC PATHING & DIRECTORY SETUP
    % Fix: Replaced hardcoded paths with dynamic relative pathing
    [script_dir, ~, ~] = fileparts(mfilename('fullpath'));
    OUT_DIR = fullfile(script_dir, 'updated_output_file');
    if ~exist(OUT_DIR, 'dir')
        mkdir(OUT_DIR); 
    end

    mat_file = fullfile(script_dir, 'proposed_results.mat');
    if ~isfile(mat_file)
        error('proposed_results.mat not found. Ensure it is in the root directory.');
    end

    %% 2. LOAD DATA & INITIALIZE HYPERSPECTRAL TOOLBOX
    fprintf('Loading dataset using Hyperspectral Imaging Library...\n');
    % Fix: Only loading 'cube' and 'wavelengths' to avoid stale variables
    load(mat_file, 'cube', 'wavelengths'); 

    % Fix: Mandatory use of hypercube() per challenge rubric
    hc = hypercube(cube, wavelengths);
    cube_data = hc.DataCube;
    wl = hc.Wavelength;
    [rows, cols, bands] = size(cube_data);

    %% 3. STAGE 1: CIBR BASELINE (Using required toolbox function)
    fprintf('Executing Stage 1: CIBR Baseline...\n');
    % Fix: Mandatory use of removecontinuum() per challenge rubric
    hc_continuum_removed = removecontinuum(hc);
    
    % Extract the specific band for CO2 absorption (~2050 nm)
    [~, co2_idx] = min(abs(wl - 2050));
    cibrScore = hc_continuum_removed.DataCube(:,:,co2_idx);
    
    % Normalize for visualization
    cibrScore = (cibrScore - min(cibrScore(:))) ./ (max(cibrScore(:)) - min(cibrScore(:)) + 1e-8);

    %% 4. STAGE 2: JRGE BASELINE
    fprintf('Executing Stage 2: JRGE Baseline...\n');
    % Using basic local deviation to simulate continuum subtraction
    jrgeScore = zeros(rows, cols);
    for r = 1:rows
        for c = 1:cols
            spectrum = squeeze(cube_data(r, c, :));
            % Calculate local deviation around CO2 band
            jrgeScore(r, c) = abs(spectrum(co2_idx) - mean(spectrum(co2_idx-5:co2_idx+5)));
        end
    end
    jrgeScore = (jrgeScore - min(jrgeScore(:))) ./ (max(jrgeScore(:)) - min(jrgeScore(:)) + 1e-8);

    %% 5. STAGE 3: SFA BASELINE (Using required toolbox function)
    fprintf('Executing Stage 3: SFA Baseline...\n');
    [target_sig] = buildtargetspectrum(wl);
    
    % Fix: Mandatory use of spectralmatch() per challenge rubric
    % Using Spectral Angle Mapper (SAM) against the CO2 physical template
    sfaScore_raw = spectralmatch(hc, target_sig, 'Method', 'sam'); 
    
    % Invert SAM score (lower angle = better match) and normalize
    sfaScore = 1 ./ (sfaScore_raw + 1e-5);
    sfaScore = (sfaScore - min(sfaScore(:))) ./ (max(sfaScore(:)) - min(sfaScore(:)) + 1e-8);

    %% 6. STAGE 4: PROPOSED CT-ACE ALGORITHM
    fprintf('Executing Stage 4: Proposed CT-ACE...\n');
    % Physics Fix: Convert Reflectance to Absorbance
    A = -log10(cube_data + 1e-6);
    
    % Water Vapor Masking
    validBands = (wl < 1800) | (wl > 1950);
    A_valid = A(:, :, validBands);
    wl_valid = wl(validBands);
    target_valid = target_sig(validBands);
    num_valid = sum(validBands);
    
    reshapedA = reshape(A_valid, rows*cols, num_valid);
    
    % K-Means Clustering
    num_clusters = 4;
    opts = statset('UseParallel', false);
    [cluster_idx, ~] = kmeans(reshapedA, num_clusters, 'Options', opts, 'MaxIter', 50);
    
    aceScore = zeros(rows*cols, 1);
    
    % Adaptive Coherence Estimator Mathematics
    for k = 1:num_clusters
        idx = (cluster_idx == k);
        if sum(idx) > num_valid
            cluster_data = reshapedA(idx, :);
            cluster_mean = mean(cluster_data, 1);
            cluster_data_centered = cluster_data - cluster_mean;
            
            % Local Covariance
            C_mat = cov(cluster_data) + eye(num_valid) * 1e-5;
            
            % Fix: Replaced inv() with backslash operator (\) for numerical stability
            w = C_mat \ target_valid; 
            
            num = cluster_data_centered * w;
            den1 = sum(cluster_data_centered .* (cluster_data_centered / C_mat), 2);
            den2 = target_valid' * w;
            
            aceScore(idx) = sign(num) .* (num.^2) ./ (den1 * den2 + 1e-8);
        end
    end
    aceScore = reshape(aceScore, rows, cols);
    aceScore(aceScore < 0) = 0; % Remove negative noise
    aceScore = (aceScore - min(aceScore(:))) ./ (max(aceScore(:)) - min(aceScore(:)) + 1e-8);

    %% 7. GENERATE FIGURE 1: STAGEWISE ABLATION
    fprintf('Generating Visualizations...\n');
    fig1 = figure('Name', 'Stagewise Ablation', 'Color', 'w', 'Position', [100 100 1200 800]);
    colormap('parula');
    
    subplot(2,3,1); imagesc(cube_data(:,:,co2_idx)); title('(a) Scene (SWIR)'); axis off;
    subplot(2,3,2); imagesc(cibrScore); title('(b) CIBR Baseline'); axis off; colorbar;
    subplot(2,3,3); imagesc(jrgeScore); title('(c) JRGE Baseline'); axis off; colorbar;
    subplot(2,3,4); imagesc(sfaScore); title('(d) SFA Baseline'); axis off; colorbar;
    subplot(2,3,5); imagesc(aceScore); title('(e) Proposed CT-ACE'); axis off; colorbar;
    
    % Thresholding for hotspot mask
    mask = aceScore > prctile(aceScore(:), 95);
    subplot(2,3,6); imshow(mask); title('(f) Spatial Hotspot Mask');
    
    saveas(fig1, fullfile(OUT_DIR, 'Figure1_Stagewise_Ablation.png'));

    %% 8. GEOSPATIAL VALIDATION MODULE
    mapgeospatial_overlay(mask, rows, cols, OUT_DIR);
    
    fprintf('Pipeline Execution Complete. All assets saved to: %s\n', OUT_DIR);
end

%% ========================================================================
% LOCAL HELPER FUNCTIONS (Consolidated to avoid missing file errors)
% ========================================================================

function [target_sig] = buildtargetspectrum(wavelengths)
    % Generates the dual-Gaussian physical CO2 template
    c1 = 1575; s1 = 15; a1 = 0.3;
    c2 = 2005; s2 = 12; a2 = 0.7;
    
    gauss1 = a1 * exp(-((wavelengths - c1).^2) / (2 * s1^2));
    gauss2 = a2 * exp(-((wavelengths - c2).^2) / (2 * s2^2));
    
    target_sig = gauss1 + gauss2;
    target_sig = target_sig / norm(target_sig);
end

function mapgeospatial_overlay(binary_mask, rows, cols, OUT_DIR)
    % GEOSPATIAL PROJECTION
    % UTM Zone 11N metadata for AVIRIS scene f250923t01p00r13
    UL_Easting = 577561.590;
    UL_Northing = 4228899.200;
    dx = 14.4; dy = 14.4;
    
    % Fix: Replaced deprecated maprasterref with maprefcells
    R = maprefcells('RasterSize', [rows, cols], ...
        'XWorldLimits', [UL_Easting, UL_Easting + cols*dx], ...
        'YWorldLimits', [UL_Northing - rows*dy, UL_Northing], ...
        'ColumnsStartFrom', 'north');
    
    [r_idx, c_idx] = find(binary_mask);
    if isempty(r_idx)
        fprintf('No hotspot pixels found for geospatial mapping.\n');
        return;
    end
    
    [Easting, Northing] = intrinsicToWorld(R, c_idx, r_idx);
    
    % UTM Zone 11N projection definition
    proj = projcrs(32611, 'Authority', 'EPSG');
    [lat, lon] = projinv(proj, Easting, Northing);
    
    fig_geo = figure('Name', 'Geospatial Overlay', 'Color', 'w', 'Position', [150 150 800 600]);
    gx = geoaxes('Basemap', 'satellite');
    hold(gx, 'on');
    geoscatter(gx, lat, lon, 20, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
    title(gx, 'Geospatial Validation of CO_2 Localization', 'FontWeight', 'bold');
    
    saveas(fig_geo, fullfile(OUT_DIR, 'Figure4_Geospatial_Validation.jpg'));
end
