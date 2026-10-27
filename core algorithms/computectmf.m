function [plume_mask, snr_map] = computectmf(radiance_cube, target_signature, background_cov, threshold)
% COMPUTECTMF Computes the Constrained Target Matching Filter (CTMF) 
% for Hyperspectral CO2 Plume Detection.
%
% Syntax:
%   [plume_mask, snr_map] = computectmf(radiance_cube, target_signature, background_cov, threshold)
%
% Inputs:
%   radiance_cube    - 3D numeric array (Height x Width x Bands) of calibrated radiance
%   target_signature - 1D vector (Bands x 1) representing the CO2 absorption spectral signature
%   background_cov   - 2D matrix (Bands x Bands) of background clutter covariance
%   threshold        - Scalar detection threshold for binary mask generation
%
% Outputs:
%   plume_mask       - 2D logical matrix (Height x Width) indicating detected CO2 enhancement
%   snr_map          - 2D floating-point map of matched filter detection scores

    % Input validation
    validateattributes(radiance_cube, {'numeric'}, {'ndim', 3}, mfilename, 'radiance_cube', 1);
    validateattributes(target_signature, {'numeric'}, {'vector'}, mfilename, 'target_signature', 2);
    validateattributes(background_cov, {'numeric'}, {'square'}, mfilename, 'background_cov', 3);
    
    [height, width, bands] = size(radiance_cube);
    
    % Ensure target vector orientation
    t = target_signature(:);
    if length(t) ~= bands
        error('Target signature length must match the number of spectral bands.');
    end
    
    % Regularize covariance matrix for numerical stability
    reg_lambda = 1e-6;
    cov_reg = background_cov + reg_lambda * eye(bands);
    
    % Pre-compute inverse covariance matrix and denominator weight
    try
        inv_cov = inv(cov_reg);
    catch
        inv_cov = pinv(cov_reg);
    end
    
    denominator = t' * inv_cov * t;
    if denominator == 0
        error('Target signature is orthogonal to inverse covariance space.');
    end
    
    % Reshape 3D radiance cube into 2D matrix (Pixels x Bands) for vectorized processing
    pixel_spectra = reshape(radiance_cube, height * width, bands);
    
    % Estimate global radiance mean for background subtraction
    mu = mean(pixel_spectra, 1);
    centered_spectra = pixel_spectra - mu;
    
    % Compute CTMF score vector across all pixels
    % CTMF = (t' * inv(K) * (x - mu)) / (t' * inv(K) * t)
    numerator_map = centered_spectra * (inv_cov * t);
    snr_vector = numerator_map / denominator;
    
    % Reshape back to 2D spatial grid
    snr_map = reshape(snr_vector, height, width);
    
    % Generate binary detection mask
    plume_mask = snr_map > threshold;
end
