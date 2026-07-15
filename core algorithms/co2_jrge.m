function [gas_map, hotspot_mask] = co2_jrge(mat_file)
% CO2_JRGE  Stage 2 — Joint Reflectance and Gas Estimator
%
% Separates CO2 gas signal from surface reflectance using an iterative
% cubic-smoothing-spline background model over the 2000–2100 nm window.
% Uses MATLAB Hyperspectral Imaging Toolbox:
%   hypercube()  — data container and wavelength accessor
%
% USAGE
%   [gas_map, hotspot_mask] = co2_jrge()
%   [gas_map, hotspot_mask] = co2_jrge(mat_file)
%
% INPUT
%   mat_file  (optional) — full path to proposed_results.mat
%
% OUTPUT
%   gas_map      [rows x cols]  Normalised gas density map in [0,1]
%   hotspot_mask [rows x cols]  Binary hotspot mask (Otsu, area >= 10 px)
%
% ALGORITHM
%   For each pixel i, a cubic smoothing spline ŝ_i^(0) = S_p[r_i, λ]
%   is fitted with smoothing parameter p = 0.90.  The gas density
%   increment at each iteration t is:
%     ĝ_i^(t) = median( ŝ_i^(t)[B] ) - median( r_i[B] )
%   where B = {b : λ_b ∈ [2000, 2100] nm}.  The spline is updated:
%     ŝ_i^(t+1) = ŝ_i^(t) + α (r_i - ŝ_i^(t)),  α = 0.05
%   over T = 3 iterations.  gasCol = Σ_t ĝ_i^(t).
%
% REFERENCES
%   Thompson et al. (2015); Marion et al. (2004)
% =========================================================================

    clc; close all;

    %% 1. RESOLVE DATA PATH
    if nargin < 1 || isempty(mat_file)
        this_dir = fileparts(mfilename('fullpath'));
        mat_file = fullfile(this_dir, '..', 'proposed_results.mat');
        mat_file = char(java.io.File(mat_file).getCanonicalPath());
    end
    if ~isfile(mat_file)
        error('co2_jrge: proposed_results.mat not found at:\n  %s', mat_file);
    end

    %% 2. LOAD DATA VIA hypercube()
    fprintf('[JRGE] Loading data...\n');
    S = load(mat_file, 'cube', 'wavelengths');
    hcube     = hypercube(double(S.cube), S.wavelengths(:));
    wavelengths = hcube.Wavelength;
    [rows, cols, ~] = size(hcube.DataCube);
    N = rows * cols;

    %% 3. BAND SELECTION: CO2 absorption window 2000–2100 nm
    co2Mask = wavelengths >= 2000 & wavelengths <= 2100;
    wlW     = wavelengths(co2Mask);
    nBw     = sum(co2Mask);
    % Normalised wavelength axis [0, 1] for spline fitting
    wlN = (wlW - wlW(1)) / (wlW(end) - wlW(1) + 1e-10);
    wlN = wlN(:)';

    %% 4. ITERATIVE SPLINE-BASED GAS ESTIMATION
    fprintf('[JRGE] Running %d-iteration spline background removal...\n', 3);
    X      = double(reshape(hcube.DataCube(:,:,co2Mask), N, nBw));
    gasCol = zeros(N, 1);

    p = 0.90;   % spline smoothing parameter
    alpha = 0.05;
    T = 3;

    for t = 1:T
        % Fit cubic smoothing spline per pixel and extract residual
        for i = 1:N
            sp     = csaps(wlN, X(i,:), p);
            spline_vals = fnval(sp, wlN);
            ghat   = median(spline_vals) - median(X(i,:));
            if ghat > 0
                gasCol(i) = gasCol(i) + ghat;
            end
        end
        % Update spectra with conservative gradient step
        slope = X(:,end) - X(:,1);
        C_lin = X(:,1) + slope .* wlN;                 % linear continuum
        residual = max(C_lin - X, 0);
        tau   = alpha * sum(residual, 2) / nBw;
        shape = residual ./ (sum(residual, 2) + 1e-10);
        X     = X + tau .* shape;
    end

    %% 5. POST-PROCESSING
    gas_map = reshape(gasCol, rows, cols);
    gas_map = (gas_map - min(gas_map(:))) / ...
              (max(gas_map(:)) - min(gas_map(:)) + 1e-10);

    %% 6. HOTSPOT MASK — Otsu threshold + area opening
    tau_otsu     = graythresh(gas_map);
    hotspot_mask = bwareaopen(gas_map > tau_otsu, 10);

    fprintf('[JRGE] Coverage: %.2f%%  |  Otsu threshold: %.4f\n', ...
            100*sum(hotspot_mask(:))/numel(hotspot_mask), tau_otsu);

    %% 7. VISUALISATION
    figure('Name','Stage 2: JRGE','Color','w','Position',[120 120 700 340]);
    subplot(1,2,1);
    imagesc(gas_map); axis image; colormap(gca,'parula'); colorbar;
    title('JRGE Gas Density Map','FontWeight','bold');

    subplot(1,2,2);
    imshow(hotspot_mask);
    title(sprintf('JRGE Hotspots (Otsu = %.3f)',tau_otsu),'FontWeight','bold');
end
