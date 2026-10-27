function [kl_divergence, metrics] = histogramcomparison(background_cube, plume_detected_cube, bin_edges)
% HISTOGRAMCOMPARISON Evaluates radiometric distribution divergence between 
% background clutter and detected CO2 plume regions.
%
% Inputs:
%   background_cube     - 3D numeric array of pure background radiance
%   plume_detected_cube - 3D numeric array of isolated plume radiance pixels
%   bin_edges           - Vector of histogram bin boundaries
%
% Outputs:
%   kl_divergence       - Scalar Kullback-Leibler divergence score
%   metrics             - Struct containing overlap coefficients and means

    if nargin < 3
        bin_edges = linspace(0, 1, 50);
    end
    
    bg_flat = background_cube(:);
    pl_flat = plume_detected_cube(:);
    
    % Compute normalized probability histograms
    p_bg = histcounts(bg_flat, bin_edges, 'Normalization', 'probability');
    p_pl = histcounts(pl_flat, bin_edges, 'Normalization', 'probability');
    
    % Smooth to prevent zero-division in KL divergence
    eps_val = 1e-10;
    p_bg = p_bg + eps_val; p_bg = p_bg / sum(p_bg);
    p_pl = p_pl + eps_val; p_pl = p_pl / sum(p_pl);
    
    % Compute Kullback-Leibler Divergence (D_KL(Plume || Background))
    kl_divergence = sum(p_pl .* log(p_pl ./ p_bg));
    
    % Additional statistical separation metrics
    metrics.overlap_coefficient = sum(min([p_bg; p_pl], [], 1));
    metrics.background_mean = mean(bg_flat);
    metrics.plume_mean = mean(pl_flat);
    
    % Visualization block for baseline reporting
    figure('Name', 'Histogram Baseline Comparison', 'Color', 'w');
    bar(bin_edges(1:end-1), [p_bg', p_pl'], 'grouped');
    xlabel('Radiance Intensity / Spectral Score');
    ylabel('Normalized Probability Density');
    legend({'Background Clutter', 'Detected Plume Region'});
    title(sprintf('Baseline Distribution Analysis (KL Divergence: %.4f)', kl_divergence));
    grid on;
end
