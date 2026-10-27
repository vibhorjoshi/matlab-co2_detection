%% Main Orchestration Script – Stage-wise CO2 Detection Pipeline
% Runs CIBR → JRGE → SFA → CT-ACE, generates the four required figures,
% and saves them to updated_output_file/.
%
% Designed to run from a clean clone with only proposed_results.mat present.

clear; clc; close all;

% -------------------- Path setup --------------------
thisFile = mfilename('fullpath');
srcDir   = fileparts(thisFile);
rootDir  = fileparts(srcDir);

addpath(fullfile(rootDir, 'core algorithms'));
addpath(srcDir);

outDir = fullfile(rootDir, 'updated_output_file');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

fprintf('[INFO] Root directory : %s\n', rootDir);
fprintf('[INFO] Output directory: %s\n', outDir);

% -------------------- Data loading --------------------
matFile = fullfile(rootDir, 'proposed_results.mat');
if exist(matFile, 'file')
    S = load(matFile, 'cube', 'wavelengths');
    cube = S.cube;
    wavelengths = S.wavelengths(:);
    fprintf('[INFO] Loaded proposed_results.mat  (%d x %d x %d)\n', size(cube));
else
    fprintf('[WARN] proposed_results.mat not found – generating synthetic demo cube.\n');
    H = 120; W = 160; B = 120;
    wavelengths = linspace(1500, 2100, B)';
    cube = 0.3 + 0.1*rand(H, W, B);
    % Inject a realistic dual-band CO2 plume
    plume = zeros(H, W);
    plume(45:75, 60:100) = 1;
    plume = imgaussfilt(plume, 3);
    for b = 1:B
        atten = 0.25*exp(-0.5*((wavelengths(b)-1575)/15)^2) + ...
                0.55*exp(-0.5*((wavelengths(b)-2005)/12)^2);
        cube(:,:,b) = cube(:,:,b) .* (1 - atten*plume);
    end
end

[H, W, ~] = size(cube);

% -------------------- Stage 1 – CIBR --------------------
fprintf('[INFO] Stage 1/4  CIBR ...\n');
[cibr_map, cibr_mask] = co2_cibr(cube, wavelengths);

% -------------------- Stage 2 – JRGE --------------------
fprintf('[INFO] Stage 2/4  JRGE ...\n');
[jrge_map, jrge_mask] = co2_jrge(cube, wavelengths);

% -------------------- Stage 3 – SFA --------------------
fprintf('[INFO] Stage 3/4  SFA ...\n');
[sfa_map, sfa_mask] = co2_sfa(cube, wavelengths);

% -------------------- Stage 4 – CT-ACE (proposed) --------------------
fprintf('[INFO] Stage 4/4  CT-ACE ...\n');
[ctace_map, ctace_mask] = computectmf(cube, wavelengths);

% -------------------- Figure 1 – Stage-wise Ablation --------------------
fprintf('[INFO] Generating Figure 1 – Stagewise Ablation\n');
fig1 = figure('Name','Stagewise Ablation','Position',[50 50 1200 700],'Color','w');

[~, iSWIR] = min(abs(wavelengths-2005));
subplot(2,3,1); imagesc(cube(:,:,iSWIR)); axis image off; colormap(gca,gray);
title('(a) Scene (SWIR)'); colorbar;

subplot(2,3,2); imagesc(cibr_map);  axis image off; colormap(gca,parula); title('(b) CIBR'); colorbar;
subplot(2,3,3); imagesc(jrge_map);  axis image off; colormap(gca,parula); title('(c) JRGE'); colorbar;
subplot(2,3,4); imagesc(sfa_map);   axis image off; colormap(gca,parula); title('(d) SFA');  colorbar;
subplot(2,3,5); imagesc(ctace_map); axis image off; colormap(gca,parula); title('(e) CT-ACE (Proposed)'); colorbar;
subplot(2,3,6); imagesc(ctace_mask); axis image off; colormap(gca,gray); title('(f) Final Hotspot Mask'); colorbar;

print(fig1, fullfile(outDir, 'Figure1_Stagewise_Ablation.png'), '-dpng', '-r150');

% -------------------- Figure 2 – Profile Analysis --------------------
fprintf('[INFO] Generating Figure 2 – Profile Analysis\n');
stats = regionprops(ctace_mask, 'Centroid');
if ~isempty(stats)
    row = round(stats(1).Centroid(2));
else
    row = round(H/2);
end
row = max(1, min(H, row));

fig2 = figure('Name','Profile Analysis','Position',[100 100 900 400],'Color','w');
plot(cibr_map(row,:),  'Color',[0.5 0.5 0.5], 'LineWidth',1.2); hold on;
plot(jrge_map(row,:),  'Color',[0.2 0.4 0.8], 'LineWidth',1.2);
plot(sfa_map(row,:),   'Color',[0.9 0.5 0.1], 'LineWidth',1.2);
plot(ctace_map(row,:), 'Color',[0.85 0.1 0.1], 'LineWidth',2.0);
xlabel(sprintf('Column index (row = %d)', row));
ylabel('Normalized detection score');
legend({'CIBR','JRGE','SFA','CT-ACE (Proposed)'}, 'Location','best');
title('Cross-section through detected plume centroid');
grid on; box on;
print(fig2, fullfile(outDir, 'Figure2_Profile_Analysis.png'), '-dpng', '-r150');

% -------------------- Figure 3 – Threshold Sensitivity --------------------
fprintf('[INFO] Generating Figure 3 – Threshold Sensitivity\n');
pcts = 85:1:99;
nPix_ctace = zeros(size(pcts));
nPix_cibr  = zeros(size(pcts));
for i = 1:numel(pcts)
    nPix_ctace(i) = sum(ctace_map(:) > prctile(ctace_map(:), pcts(i)));
    nPix_cibr(i)  = sum(cibr_map(:)  > prctile(cibr_map(:),  pcts(i)));
end
fig3 = figure('Name','Threshold Sensitivity','Position',[150 150 700 400],'Color','w');
plot(pcts, nPix_ctace, '-o', 'LineWidth',1.5, 'MarkerFaceColor',[0.85 0.1 0.1]); hold on;
plot(pcts, nPix_cibr,  '-s', 'LineWidth',1.5, 'MarkerFaceColor',[0.5 0.5 0.5]);
xlabel('Percentile threshold');
ylabel('Number of pixels above threshold');
legend({'CT-ACE','CIBR'}, 'Location','best');
title('Threshold sensitivity');
grid on; box on;
print(fig3, fullfile(outDir, 'Figure3_Threshold_Sensitivity.png'), '-dpng', '-r150');

% -------------------- Figure 4 – Geospatial Validation --------------------
fprintf('[INFO] Generating Figure 4 – Geospatial Validation\n');
fig4 = figure('Name','Geospatial Validation','Position',[200 200 700 500],'Color','w');
% Approximate UTM Zone 11N coordinates (matches README / original flightline)
UL_E = 577561.59; UL_N = 4228899.2; dx = 14.4; dy = 14.4;
x = UL_E + (0:W-1)*dx;
y = UL_N - (0:H-1)*dy;
imagesc(x, y, ctace_mask);
axis image xy;
colormap(gray);
xlabel('Easting (m, UTM Zone 11N)');
ylabel('Northing (m, UTM Zone 11N)');
title('CT-ACE Hotspot Mask – Geospatial Overlay');
colorbar;
print(fig4, fullfile(outDir, 'Figure4_Geospatial_Validation.png'), '-dpng', '-r150');

% -------------------- Baseline histogram comparison --------------------
fprintf('[INFO] Running histogram comparison (baseline vs CT-ACE)\n');
try
    histogram_comparison(cibr_map, ctace_map, outDir);
catch ME
    warning('histogram_comparison failed: %s', ME.message);
end

fprintf('\n========================================\n');
fprintf(' PIPELINE COMPLETED SUCCESSFULLY\n');
fprintf(' CIBR  hotspots : %d pixels\n', sum(cibr_mask(:)));
fprintf(' JRGE  hotspots : %d pixels\n', sum(jrge_mask(:)));
fprintf(' SFA   hotspots : %d pixels\n', sum(sfa_mask(:)));
fprintf(' CT-ACE hotspots: %d pixels\n', sum(ctace_mask(:)));
fprintf(' Figures saved to: %s\n', outDir);
fprintf('========================================\n');
