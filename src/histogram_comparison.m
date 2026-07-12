% =========================================================================
% histogram_comparison.m
% STEP 4 — Score Distribution Comparison (Computed from Scratch)
%
% Proves background clutter suppression by comparing:
%   • Baseline: Standard Global Matched Filter (Raw Reflectance)
%   • Proposed: CT-ACE (Absorbance Preprocessing + Clustering)
% =========================================================================
clear; clc; close all;

fprintf('\n%s\n',  repmat('=',1,60));
fprintf(' STEP 4  |  SCORE HISTOGRAM COMPARISON\n');
fprintf('%s\n\n', repmat('=',1,60));

%% 1. LOAD RAW DATA
mat_file = fullfile(pwd, 'proposed_results.mat');
if ~isfile(mat_file)
    error('proposed_results.mat not found. Ensure it is in your repository.');
end
load(mat_file, 'cube', 'wavelengths'); 
wavelengths = wavelengths(:);
[rows, cols, ~] = size(cube);

% Valid SWIR Band Selection
valid_bands = (wavelengths >= 1500 & wavelengths <= 2100) & ~(wavelengths >= 1800 & wavelengths <= 1950);
num_valid = sum(valid_bands);
wl_valid = wavelengths(valid_bands);

% Target Signature Generation
amp1 = 0.30; cen1 = 1575; sig1 = 15;
amp2 = 0.70; cen2 = 2005; sig2 = 12;
co2_sig = zeros(num_valid, 1);
for b = 1:num_valid
    lam = wl_valid(b);
    co2_sig(b) = amp1*exp(-0.5*((lam-cen1)/sig1)^2) + amp2*exp(-0.5*((lam-cen2)/sig2)^2);
end
co2_sig = co2_sig / norm(co2_sig);

%% 2. COMPUTE BASELINE: Global Matched Filter (Raw Reflectance)
fprintf('Computing Baseline Matched Filter...\n');
X_raw = double(reshape(cube(:,:,valid_bands), rows*cols, num_valid));
mu_raw = mean(X_raw, 1)';
X_c_raw = X_raw - mu_raw';

C_global = cov(X_raw) + eye(num_valid)*1e-5;
invC_global = inv(C_global);
w_base = invC_global * co2_sig;

bs_raw = X_c_raw * w_base;
bs_raw(bs_raw < 0) = 0; % Rectify
bs = bs_raw / (max(bs_raw(:)) + 1e-10); % Normalize [0,1]

%% 3. COMPUTE PROPOSED: CT-ACE (Absorbance + Clustering)
fprintf('Computing Proposed CT-ACE...\n');
A_cube = -log(max(double(cube(:,:,valid_bands)), 1e-5)); % Absorbance Transform
reshapedA = reshape(A_cube, rows*cols, num_valid);

rng(42);
[clusterIdx, ~] = kmeans(reshapedA, 4, 'MaxIter', 100, 'Distance', 'sqeuclidean');
cs_raw = zeros(rows*cols, 1);

for k = 1:4
    idx = (clusterIdx == k);
    X_c = reshapedA(idx, :) - mean(reshapedA(idx, :), 1);
    if sum(idx) < num_valid + 5, continue; end
    
    invC = inv(cov(reshapedA(idx, :)) + eye(num_valid)*1e-5);
    w = invC * co2_sig;
    
    num = X_c * w; 
    den1 = sum(X_c .* (X_c * invC), 2); 
    den2 = co2_sig' * w; 
    cs_raw(idx) = sign(num) .* (num.^2) ./ (den1 * den2 + 1e-8);
end

cs_raw(cs_raw < 0) = 0;
cs_raw = medfilt2(reshape(cs_raw, rows, cols), [3 3]); % Spatial filter
cs = cs_raw(:) / (max(cs_raw(:)) + 1e-10); % Normalize [0,1]

%% 4. DESCRIPTIVE STATISTICS
stat_names = {'Mean', 'Median', 'StdDev', 'Maximum', '95th Pctile'};
st_b = [mean(bs), median(bs), std(bs), max(bs), prctile(bs, 95)];
st_c = [mean(cs), median(cs), std(cs), max(cs), prctile(cs, 95)];

fprintf('\n%-18s  %-14s  %-14s\n', 'Metric', 'Baseline CTMF', 'Proposed CT-ACE');
fprintf('%s\n', repmat('-',50,1));
for k = 1 : 5
    fprintf('%-18s  %-14.5f  %-14.5f\n', stat_names{k}, st_b(k), st_c(k));
end
fprintf('\n');

%% 5. FIGURE GENERATION
all_vals  = [bs; cs];
lo_edge   = prctile(all_vals, 0.5);
hi_edge   = prctile(all_vals, 99.5);
N_BINS    = 120;
edges     = linspace(lo_edge, hi_edge, N_BINS + 1);

fig = figure('Name','Score Histograms', 'Color', 'w', 'Position',[100 100 800 500]);
hold on;

h1 = histogram(bs, edges, ...
    'Normalization','probability', ...
    'FaceColor',  [0.22 0.47 0.74], ...
    'FaceAlpha',  0.65, ...
    'EdgeColor',  'none', ...
    'DisplayName','Baseline (Raw Reflectance)');

h2 = histogram(cs, edges, ...
    'Normalization','probability', ...
    'FaceColor',  [0.84 0.19 0.13], ...
    'FaceAlpha',  0.65, ...
    'EdgeColor',  'none', ...
    'DisplayName','Proposed CT-ACE (Conditioned)');

% Mean indicator lines
xline(mean(bs), '--', 'Color', [0.22 0.47 0.74], 'LineWidth', 2, ...
    'Label', sprintf('\\mu_{base}=%.3f', mean(bs)), ...
    'LabelVerticalAlignment','bottom', 'HandleVisibility','off');
xline(mean(cs), '--', 'Color', [0.84 0.19 0.13], 'LineWidth', 2, ...
    'Label', sprintf('\\mu_{prop}=%.3f', mean(cs)), ...
    'LabelVerticalAlignment','top', 'HandleVisibility','off');

hold off;
legend([h1, h2], 'Location','northeast', 'FontSize', 10);
xlabel('Normalized Anomaly Score', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Probability', 'FontSize', 12, 'FontWeight', 'bold');
title('Score Distribution: Background Clutter Suppression', 'FontSize', 14, 'FontWeight', 'bold');

% Statistics annotation box
ann_str = sprintf('Baseline  —  Max = %.3f   95th = %.3f   σ = %.4f\n', st_b(4), st_b(5), st_b(3));
ann_str = [ann_str, sprintf('Proposed  —  Max = %.3f   95th = %.3f   σ = %.4f', st_c(4), st_c(5), st_c(3))];
annotation('textbox', [0.15 0.70 0.40 0.12], ...
    'String',          ann_str, ...
    'FitBoxToText',    'off', ...
    'BackgroundColor', [1 1 0.96], ...
    'EdgeColor',       [0.6 0.6 0.6], ...
    'FontSize',        9);

grid on; box on;

OUT_DIR = 'output_figures'; if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end
OUT_PNG = fullfile(OUT_DIR, 'histogram_scores.png');
saveas(fig, OUT_PNG);
fprintf('Saved → %s\n', OUT_PNG);