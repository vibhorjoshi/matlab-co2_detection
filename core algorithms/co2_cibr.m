function [cibr_map, plume_mask] = co2cibr(radiance_cube, wavelengths, absorption_band_idx, ref_band_indices, threshold)
% CO2CIBR Computes the Carbon Dioxide Absorption Band Index Ratio (CIBR).
%
% Inputs:
%   radiance_cube        - 3D numeric array (Height x Width x Bands)
%   wavelengths          - 1D vector of spectral band centers
%   absorption_band_idx  - Index of the primary CO2 absorption feature (e.g., ~2.0 µm region)
%   ref_band_indices     - 1 or 2 indices representing reference continuum bands
%   threshold            - Scalar threshold for plume delineation
%
% Outputs:
%   cibr_map             - 2D spatial ratio map
%   plume_mask           - 2D logical mask

    validateattributes(radiance_cube, {'numeric'}, {'ndim', 3}, mfilename, 'radiance_cube', 1);
    [height, width, bands] = size(radiance_cube);
    
    pixel_spectra = reshape(radiance_cube, height * width, bands);
    
    % Extract absorption and reference radiances
    r_abs = pixel_spectra(:, absorption_band_idx);
    
    if length(ref_band_indices) == 2
        r_ref = mean(pixel_spectra(:, ref_band_indices), 2);
    else
        r_ref = pixel_spectra(:, ref_band_indices);
    end
    
    % CIBR Formulation: Ratio of absorption band to interpolated continuum reference
    ratio_vector = r_abs ./ max(r_ref, 1e-4);
    
    cibr_map = reshape(ratio_vector, height, width);
    
    % Low ratio values indicate enhanced absorption (CO2 presence)
    plume_mask = cibr_map < threshold;
end
