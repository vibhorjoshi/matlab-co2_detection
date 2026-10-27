function [sfa_map, plume_mask] = co2sfa(radiance_cube, threshold)
% CO2SFA Applies Slow Feature Analysis to isolate invariant/slowly varying components.
%
% Inputs:
%   radiance_cube - 3D numeric array (Height x Width x Bands)
%   threshold     - Statistical anomaly threshold
%
% Outputs:
%   sfa_map       - 2D slow feature projection map
%   plume_mask    - 2D binary plume mask

    [height, width, bands] = size(radiance_cube);
    X = double(reshape(radiance_cube, height * width, bands));
    
    % Zero-mean normalization
    X_mean = mean(X, 1);
    X_centered = X - X_mean;
    
    % Covariance matrix and whitening
    C = (X_centered' * X_centered) / (size(X, 1) - 1);
    [V, D] = eig(C);
    
    evals = diag(D);
    [evals_sorted, sort_idx] = sort(evals, 'descend');
    V_sorted = V(:, sort_idx);
    
    % Whitening matrix
    W = diag(1 ./ sqrt(evals_sorted + 1e-5)) * V_sorted';
    X_white = X_centered * W';
    
    % Time derivative approximation (spatial finite differences)
    X_spatial = reshape(X_white, height, width, bands);
    dX_dx = diff(X_spatial, 1, 2); % horizontal derivative
    dX_dy = diff(X_spatial, 1, 1); % vertical derivative
    
    % SFA optimization via derivative covariance eigen-decomposition
    dX_flat = [reshape(dX_dx, [], bands); reshape(dX_dy(1:end-1,:,:), [], bands)];
    Omega = (dX_flat' * dX_flat) / size(dX_flat, 1);
    
    [W_sfa, ~] = eig(Omega);
    sfa_projection = X_white * W_sfa(:, 1); % Slowest feature
    
    sfa_map = reshape(sfa_projection, height, width);
    
    % Anomaly thresholding based on robust median absolute deviation (MAD)
    med_val = median(sfa_map(:));
    mad_val = mad(sfa_map(:), 1);
    plume_mask = abs(sfa_map - med_val) > (threshold * mad_val);
end
