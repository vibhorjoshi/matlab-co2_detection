function [jrge_score_map, plume_mask] = co2jrge(radiance_cube, background_subspace, target_vector, threshold)
% CO2JRGE Computes Joint Residual Gas Enhancement detection scores.
%
% Inputs:
%   radiance_cube       - 3D numeric array (Height x Width x Bands)
%   background_subspace - Orthonormal basis matrix of background clutter (Bands x SubspaceDim)
%   target_vector       - Target gas absorption vector (Bands x 1)
%   threshold           - Detection threshold
%
% Outputs:
%   jrge_score_map      - 2D spatial score map
%   plume_mask          - 2D logical detection mask

    [height, width, bands] = size(radiance_cube);
    pixel_spectra = reshape(radiance_cube, height * width, bands);
    
    t = target_vector(:);
    
    % Project spectra onto background orthogonal complement space
    P_orth = eye(bands) - (background_subspace * background_subspace');
    
    residual_spectra = pixel_spectra * P_orth;
    projected_target = P_orth * t;
    
    norm_proj_target = norm(projected_target);
    if norm_proj_target > 0
        unit_target = projected_target / norm_proj_target;
    else
        unit_target = projected_target;
    end
    
    score_vector = residual_spectra * unit_target;
    jrge_score_map = reshape(score_vector, height, width);
    
    plume_mask = jrge_score_map > threshold;
end
