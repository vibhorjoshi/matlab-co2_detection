function [kl_div, metrics] = histogram_comparison(baseline_map, proposed_map, outDir, labels)
%HISTOGRAM_COMPARISON  Score-distribution comparison (Baseline vs CT-ACE)
%
%   [kl_div, metrics] = histogram_comparison()                          % runs both methods
%   [kl_div, metrics] = histogram_comparison(baseline_map, proposed_map)
%   [kl_div, metrics] = histogram_comparison(baseline_map, proposed_map, outDir)
%   [kl_div, metrics] = histogram_comparison(..., outDir, {'Baseline','Proposed'})
%
% Computes KL divergence and overlap between two detection-score maps in [0,1].
% With no arguments it scores the fallback data with the global matched-filter
% BASELINE (reflectance domain, correctly signed target) and with CT-ACE, using the
% same post-processing for both.
%
% The figure is saved to <outDir>/histogram_scores.png; the default outDir is
% <repo>/updated_output_file, resolved from this file's location (never from pwd),
% so run('src/histogram_comparison.m') works from any folder.
%
% NOTE: a histogram shows that two scores are distributed differently, not that one
% is more accurate. Accuracy is measured against ground truth in main_pipeline
% (Table1 / Table2, Figure 5).

    U = co2_util();                                  % needs "core algorithms" on the path
    if nargin < 3 || isempty(outDir)
        outDir = fullfile(U.root(), 'updated_output_file');
    end
    if nargin < 4 || isempty(labels)
        labels = {'Baseline (global matched filter)', 'Proposed (CT-ACE)'};
    end
    if nargin < 2
        baseline_map = U.baselineMF();
        proposed_map = computectmf();
    end

    b = baseline_map(:);   b = b(isfinite(b));
    p = proposed_map(:);   p = p(isfinite(p));

    nBins = 50;
    edges = linspace(0, 1, nBins + 1);
    pb = bin_probability(b, nBins);
    pp = bin_probability(p, nBins);

    eps_val = 1e-10;
    pb = pb + eps_val;  pb = pb / sum(pb);
    pp = pp + eps_val;  pp = pp / sum(pp);

    kl_div = sum(pp .* log(pp ./ pb));              % D_KL(proposed || baseline)

    metrics.overlap_coefficient = sum(min(pb, pp));
    metrics.baseline_mean  = mean(b);
    metrics.proposed_mean  = mean(p);
    metrics.baseline_std   = std(b);
    metrics.proposed_std   = std(p);

    fprintf('%-12s  %-12s  %-12s\n', 'Metric', 'Baseline', 'Proposed');
    fprintf('%-12s  %-12.5f  %-12.5f\n', 'Mean',   metrics.baseline_mean, metrics.proposed_mean);
    fprintf('%-12s  %-12.5f  %-12.5f\n', 'Std dev', metrics.baseline_std, metrics.proposed_std);
    fprintf('KL(proposed||baseline) = %.4f | overlap = %.4f\n', kl_div, metrics.overlap_coefficient);

    % Visualisation
    fig = figure('Name', 'Histogram Comparison - Baseline vs CT-ACE', 'Color', 'w', 'Visible', 'off');
    bar(edges(1:end-1), [pb; pp]', 'grouped');
    set(gca, 'YScale', 'log');
    xlabel('Normalized detection score');
    ylabel('Probability (log scale)');
    legend(labels, 'Location', 'best');
    title(sprintf('Score distribution comparison  (KL = %.4f)', kl_div));
    grid on; box on;

    U.saveFig(fig, fullfile(outDir, 'histogram_scores.png'));
    close(fig);
    fprintf('Saved -> %s\n', fullfile(outDir, 'histogram_scores.png'));
end

function pr = bin_probability(x, nBins)
%BIN_PROBABILITY  Histogram probabilities on nBins equal bins over [0,1] (no toolbox needed).
    idx = min(max(floor(x * nBins) + 1, 1), nBins);
    pr  = accumarray(idx(:), 1, [nBins 1])' / numel(x);
end
