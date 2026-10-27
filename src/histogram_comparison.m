function [kl_div, metrics] = histogram_comparison(baseline_map, proposed_map, outDir)
%HISTOGRAM_COMPARISON  Score-distribution comparison (Baseline vs CT-ACE)
%
%   [kl_div, metrics] = histogram_comparison(baseline_map, proposed_map)
%   [kl_div, metrics] = histogram_comparison(baseline_map, proposed_map, outDir)
%
% Computes KL divergence and overlap between the two detection-score maps.
% When outDir is supplied the figure is also saved there.

    if nargin < 3
        outDir = [];
    end

    b = baseline_map(:);
    p = proposed_map(:);
    b = b(isfinite(b));
    p = p(isfinite(p));

    bin_edges = linspace(0, 1, 51);
    pb = histcounts(b, bin_edges, 'Normalization', 'probability');
    pp = histcounts(p, bin_edges, 'Normalization', 'probability');

    eps_val = 1e-10;
    pb = pb + eps_val; pb = pb / sum(pb);
    pp = pp + eps_val; pp = pp / sum(pp);

    kl_div = sum(pp .* log(pp ./ pb));          % D_KL(proposed || baseline)

    metrics.overlap_coefficient = sum(min(pb, pp));
    metrics.baseline_mean  = mean(b);
    metrics.proposed_mean  = mean(p);
    metrics.baseline_std   = std(b);
    metrics.proposed_std   = std(p);

    % Visualisation
    fig = figure('Name','Histogram Comparison – Baseline vs CT-ACE','Color','w');
    bar(bin_edges(1:end-1), [pb; pp]', 'grouped');
    xlabel('Normalized detection score');
    ylabel('Probability density');
    legend({'Baseline (CIBR)','Proposed (CT-ACE)'}, 'Location','best');
    title(sprintf('Score distribution comparison  (KL = %.4f)', kl_div));
    grid on; box on;

    if ~isempty(outDir) && exist(outDir, 'dir')
        print(fig, fullfile(outDir, 'histogram_scores.png'), '-dpng', '-r150');
    end
end
