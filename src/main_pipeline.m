%% Main Orchestration Script for Hyperspectral CO2 Detection Pipeline
% Automatically manages workspace paths, executes core detectors, 
% and runs comparative baseline validation.

clear; clc; close all;

% --- 1. Automated Path Initialization ---
current_file_path = mfilename('fullpath');
[src_dir, ~, ~] = fileparts(current_file_path);
root_dir = fileparts(src_dir);

addpath(fullfile(root_dir, 'core'));
addpath(fullfile(root_dir, 'src'));

fprintf('[INFO] Root directory initialized: %s\n', root_dir);

% --- 2. Synthetic/Sample Hyperspectral Data Generation (Mock HSI Cube) ---
% In production, replace with ENVI .hdr/.bsq loaders (e.g., multibandread)
height = 100; width = 100; bands = 150;
wavelengths = linspace(2.0, 2.5, bands); % Simulated shortwave infrared (SWIR) region

fprintf('[INFO] Generating synthetic radiance hyperspectral cube (%dx%dx%d)...\n', height, width, bands);
radiance_cube = rand(height, width, bands) * 0.5 + 0.1;

% Inject a synthetic CO2 plume enhancement signature in center pixels
synthetic_plume_sig = exp(-((wavelengths - 2.006).^2) / (2 * 0.01^2))'; 
for r = 40:60
    for c = 40:60
        radiance_cube(r, c, :) = squeeze(radiance_cube(r, c, :)) + 0.3 * synthetic_plume_sig;
    end
end

% --- 3. Execute Core Detectors ---
fprintf('[INFO] Executing core detection algorithms...\n');

% A. Constrained Target Matching Filter (CTMF)
background_cov = cov(reshape(radiance_cube, height * width, bands));
[ctmf_mask, ctmf_scores] = computectmf(radiance_cube, synthetic_plume_sig, background_cov, 2.5);

% B. Carbon Dioxide Absorption Band Index Ratio (CIBR)
[cibr_map, cibr_mask] = co2cibr(radiance_cube, wavelengths, 50, [45, 55], 0.85);

% C. Slow Feature Analysis (SFA)
[sfa_map, sfa_mask] = co2sfa(radiance_cube, 3.0);

% --- 4. Baseline Validation via Histogram Comparison ---
fprintf('[INFO] Running statistical baseline validation...\n');
bg_sample = radiance_cube(1:20, 1:20, 50);
plume_sample = radiance_cube(45:55, 45:55, 50);
[kl_div, baseline_metrics] = histogramcomparison(bg_sample, plume_sample);

fprintf('========================================\n');
fprintf(' PIPELINE EXECUTION COMPLETED SUCCESSFULLY\n');
fprintf(' KL Divergence Score : %.4f\n', kl_div);
fprintf(' CTMF Detected Pixels: %d\n', sum(ctmf_mask(:)));
fprintf(' CIBR Detected Pixels: %d\n', sum(cibr_mask(:)));
fprintf('========================================\n');

% --- 5. Visualization Summary ---
figure('Name', 'CO2 Plume Detection Results', 'Position', [100, 100, 1200, 400]);

subplot(1, 3, 1);
imagesc(ctmf_scores);
colorbar;
title('CTMF Detection Score Map');
xlabel('X (pixels)'); ylabel('Y (pixels)');

subplot(1, 3, 2);
imagesc(cibr_map);
colorbar;
title('CIBR Ratio Map');
xlabel('X (pixels)'); ylabel('Y (pixels)');

subplot(1, 3, 3);
imagesc(ctmf_mask | cibr_mask);
colormap(gca, gray);
title('Fused Plume Detection Mask');
xlabel('X (pixels)'); ylabel('Y (pixels)');
