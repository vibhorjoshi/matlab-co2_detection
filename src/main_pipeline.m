function results = main_pipeline(varargin)
%MAIN_PIPELINE  Stage-wise CO2 detection pipeline: CIBR -> JRGE -> SFA -> CT-ACE (+ baseline).
%
%   results = main_pipeline()                          fallback data (or synthetic scene)
%   results = main_pipeline(cube, wavelengths)         in-memory cube
%   results = main_pipeline('file.mat')
%   results = main_pipeline('scene.img', 'scene.hdr', 'Crop', [r1 r2 c1 c2])
%   results = main_pipeline(..., 'OutDir', folder, 'Validate', false, 'Seeds', 1:5)
%
% The pipeline CALLS co2_cibr, co2_jrge, co2_sfa and computectmf (and the global
% matched-filter baseline) - nothing is re-implemented here. Every path is resolved
% from this file's location, never from pwd, so it runs from any folder:
%     >> run('src/main_pipeline.m')      or      >> main_pipeline
%
% Outputs (default folder <repo>/updated_output_file):
%   Figure1_Stagewise_Ablation.png      Figure2_Profile_Analysis.png
%   Figure3_Threshold_Sensitivity.png   Figure4_Geospatial_Validation.png (real scenes)
%   Table1_Selectivity_Metrics.csv      histogram_scores.png
%   Table2_Synthetic_Validation.csv     Figure5_Synthetic_Validation.png  ('Validate', true)
 
    % -------------------- Path setup (location-independent) --------------------
    rootDir = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(rootDir, 'core algorithms'));
    addpath(fullfile(rootDir, 'src'));
    U = co2_util();
 
    % -------------------- Options --------------------
    outDir = fullfile(rootDir, 'updated_output_file');
    doValidate = true;  seeds = 1:3;  amps = [0.01 0.02 0.04 0.08];
    rest = {};  i = 1;
    while i <= numel(varargin)
        a = varargin{i};
        if ischar(a) && i < numel(varargin) && any(strcmpi(a, {'OutDir', 'Validate', 'Seeds', 'Amplitudes'}))
            switch lower(a)
                case 'outdir',     outDir = varargin{i+1};
                case 'validate',   doValidate = logical(varargin{i+1});
                case 'seeds',      seeds = varargin{i+1};
                case 'amplitudes', amps = varargin{i+1};
            end
            i = i + 2;
        else
            rest{end+1} = a;  i = i + 1; %#ok<AGROW>
        end
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    fprintf('[INFO] Root directory : %s\n', rootDir);
    fprintf('[INFO] Output directory: %s\n', outDir);
 
    % -------------------- Data loading (validated, header-based georeference) -----
    [D, ~, extra] = U.prepare(rest);
    if ~isempty(extra), error('main_pipeline:args', 'Unrecognised arguments.'); end
    fprintf('[INFO] Source: %s  (%d x %d x %d) | %s\n', D.source, size(D.cube), D.report.message);
    if ~isempty(D.truth), fprintf('[INFO] Scene carries injected ground truth.\n'); end
 
    % -------------------- Stages --------------------
    t0 = tic;
    fprintf('[INFO] Stage 1/4  CIBR ...\n');   [s, m, i1] = co2_cibr(D);        R.cibr  = struct('score', s, 'mask', m, 'info', i1);
    fprintf('[INFO] Stage 2/4  JRGE ...\n');   [s, m, i2] = co2_jrge(D);        R.jrge  = struct('score', s, 'mask', m, 'info', i2);
    fprintf('[INFO] Stage 3/4  SFA ...\n');    [s, m, i3] = co2_sfa(D);         R.sfa   = struct('score', s, 'mask', m, 'info', i3);
    fprintf('[INFO] Stage 4/4  CT-ACE ...\n'); [s, m, i4] = computectmf(D);     R.ctace = struct('score', s, 'mask', m, 'info', i4);
    fprintf('[INFO] Baseline   global MF ...\n'); [s, m, i5] = U.baselineMF(D); R.base  = struct('score', s, 'mask', m, 'info', i5);
    fprintf('[INFO] Detectors finished in %.1f s\n', toc(t0));
 
    % -------------------- Figures 1-3 and Table 1 --------------------
    make_figures(R, D, outDir, U);
 
    % -------------------- Figure 4: geospatial overlay (header georeference) -----
    U.geoOverlay(R.ctace.mask, D.georef, fullfile(outDir, 'Figure4_Geospatial_Validation.png'));
 
    % -------------------- Score histograms --------------------
    histogram_comparison(R.base.score, R.ctace.score, outDir);
 
    % -------------------- Ground-truth validation (synthetic injection) ----------
    if doValidate
        run_validation_study(U, outDir, seeds, amps);
    end
    fprintf('[INFO] Done. Outputs in %s\n', outDir);
    results = R;
end
 
% =========================================================================
function make_figures(R, D, outDir, U)
%MAKE_FIGURES  Generates Figures 1-3 and Table 1 from the detector outputs.
%
%   R  struct with fields cibr, jrge, sfa, ctace, base; each a struct with
%      .score .mask .info (as returned by the detectors)
%   D  data struct (for the scene panel and, if present, ground truth)
%
%   Every committed figure/table is produced here:
%       Figure1_Stagewise_Ablation.png    Figure2_Profile_Analysis.png
%       Figure3_Threshold_Sensitivity.png Table1_Selectivity_Metrics.csv
    if exist(outDir, 'dir') ~= 7, mkdir(outDir); end
    names = {'cibr', 'jrge', 'sfa', 'ctace', 'base'};
    labels = {'CIBR', 'JRGE', 'SFA', 'CT-ACE (proposed)', 'Global MF (baseline)'};
    cfg = U.cfg();
    [~, ib] = min(abs(D.wl - cfg.cibr.center));
 
    % ---- Figure 1: stagewise comparison ----------------------------------
    f1 = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1200 800]);
    subplot(2, 3, 1);  imagesc(D.cube(:, :, ib));  axis image off;  colorbar;
    title(sprintf('(a) Scene, %.0f nm', D.wl(ib)));
    for k = 1:4
        subplot(2, 3, k + 1);  imagesc(R.(names{k}).score, [0 1]);  axis image off;  colorbar;
        title(sprintf('(%c) %s', 'a' + k, labels{k}));
    end
    subplot(2, 3, 6);  imshow(R.ctace.mask);  title('(f) CT-ACE hotspot mask');
    colormap(f1, 'parula');
    U.saveFig(f1, fullfile(outDir, 'Figure1_Stagewise_Ablation.png'));  close(f1);
 
    % ---- Figure 2: profile through the largest CT-ACE hotspot ------------
    mk = R.ctace.mask;
    cc = regionprops(mk, 'Area', 'Centroid');
    if isempty(cc)
        [~, i] = max(R.ctace.score(:));  [row, ~] = ind2sub(size(mk), i);
    else
        [~, j] = max([cc.Area]);  row = round(cc(j).Centroid(2));
    end
    f2 = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1000 420]);
    hold on;  grid on;  box on;
    for k = 1:5
        plot(R.(names{k}).score(row, :), 'LineWidth', 1.3, 'DisplayName', labels{k});
    end
    if ~isempty(D.truth)
        area(double(D.truth(row, :)) * 1.0, 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'DisplayName', 'Injected plume (truth)');
    end
    xlabel('Column (pixel)');  ylabel('Normalised score');
    title(sprintf('Cross-section through row %d (centroid of largest CT-ACE hotspot)', row));
    legend('Location', 'best');
    U.saveFig(f2, fullfile(outDir, 'Figure2_Profile_Analysis.png'));  close(f2);
 
    % ---- Figure 3: threshold sensitivity ----------------------------------
    pct = 85:1:99.5;
    cov_ = zeros(5, numel(pct));  f1s = nan(5, numel(pct));
    for k = 1:5
        raw = R.(names{k}).info.raw;  ok = ~D.invalid;
        for p = 1:numel(pct)
            m = bwareaopen(raw > prctile(raw(ok), pct(p)) & ok, cfg.threshold.minArea);
            cov_(k, p) = 100 * nnz(m) / nnz(ok);
            if ~isempty(D.truth), f1s(k, p) = 2*nnz(m & D.truth) / max(nnz(m) + nnz(D.truth), 1); end
        end
    end
    f3 = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1100 420]);
    subplot(1, 2, 1);  plot(pct, cov_', 'LineWidth', 1.3);  grid on;
    xlabel('Threshold percentile');  ylabel('Hotspot coverage (%)');  legend(labels, 'Location', 'northeast');
    title('Coverage vs threshold');
    subplot(1, 2, 2);
    if ~isempty(D.truth)
        plot(pct, f1s', 'LineWidth', 1.3);  ylabel('F1 against injected plume');
        title('F1 vs threshold (ground truth)');
    else
        text(0.5, 0.5, 'No ground truth for this scene: F1 not defined', 'HorizontalAlignment', 'center');
        axis off;
    end
    grid on;  xlabel('Threshold percentile');
    U.saveFig(f3, fullfile(outDir, 'Figure3_Threshold_Sensitivity.png'));  close(f3);
 
    % ---- Table 1 -----------------------------------------------------------
    fid = fopen(fullfile(outDir, 'Table1_Selectivity_Metrics.csv'), 'w');
    fprintf(fid, 'Algorithm,Coverage_pct,Mean_Score_In_Mask,Mean_Score_Outside,Self_SBR,AUC_vs_truth,F1_vs_truth,Reference\n');
    for k = 1:5
        s = R.(names{k}).score;  m = R.(names{k}).mask;
        sbr = mean(s(m)) / (mean(s(~m)) + eps);
        if ~isempty(D.truth)
            mt = U.metrics(R.(names{k}).info.raw, D.truth, m);
            auc = sprintf('%.4f', mt.auc);  f1 = sprintf('%.4f', mt.f1);
            ref = 'injected-plume ground truth (AUC/F1); SBR is self-referential';
        else
            auc = 'NA';  f1 = 'NA';
            ref = 'no ground truth: Self_SBR is DESCRIPTIVE ONLY (mask defines target and background)';
        end
        fprintf(fid, '"%s",%.2f,%.4f,%.4f,%.3f,%s,%s,"%s"\n', labels{k}, ...
                R.(names{k}).info.coverage, mean(s(m)), mean(s(~m)), sbr, auc, f1, ref);
    end
    fclose(fid);
end
 
% =========================================================================
function T = run_validation_study(U, outDir, seeds, amps)
%RUN_VALIDATION_STUDY  Ground-truth validation on scenes with a KNOWN injected plume.
%
%   Every method is scored on the SAME synthetic scenes (U.synth): threshold-free AUC
%   (tie-corrected) and probability of detection at Pfa = 1e-2, averaged over seeds.
%       Global MF (baseline)      reflectance, one covariance
%       ACE, absorbance, K=1      ablation: clustering removed
%       ACE, reflectance, K=4     ablation: absorbance removed
%       CT-ACE (proposed)         absorbance + K=4 clusters
%       CIBR / JRGE / SFA         Stages 1-3
%   Writes Table2_Synthetic_Validation.csv and Figure5_Synthetic_Validation.png.
%   Scope: a controlled test bench; it does not replace validation on a real
%   flightline over a known emitter (see README, "Validation").
    cK1  = struct('ctace', struct('K', 1));
    cRef = struct('ctace', struct('absorbance', false));
    methods = { ...
        'Global MF (baseline)',    @(D) U.baselineMF(D);
        'ACE, absorbance, K=1',    @(D) computectmf(D, 'Config', cK1);
        'ACE, reflectance, K=4',   @(D) computectmf(D, 'Config', cRef);
        'CT-ACE (proposed)',       @(D) computectmf(D);
        'CIBR',                    @(D) co2_cibr(D);
        'JRGE',                    @(D) co2_jrge(D);
        'SFA',                     @(D) co2_sfa(D)};
    nM = size(methods, 1);  nA = numel(amps);  nS = numel(seeds);
    AUC = zeros(nM, nA, nS);  PD = zeros(nM, nA, nS);
    for a = 1:nA
        for s = 1:nS
            D = U.synth('Amplitude', amps(a), 'Seed', seeds(s));
            for m = 1:nM
                [~, ~, info] = methods{m, 2}(D);
                mt = U.metrics(info.raw, D.truth);
                AUC(m, a, s) = mt.auc;  PD(m, a, s) = mt.pd_pfa_1e2;
            end
        end
        fprintf('[validation] amplitude %.3f done\n', amps(a));
    end
 
    T = struct('method', {methods(:, 1)}, 'amplitude', amps, ...
               'aucMean', mean(AUC, 3), 'aucStd', std(AUC, 0, 3), 'pdMean', mean(PD, 3));
    if exist(outDir, 'dir') ~= 7, mkdir(outDir); end
    fid = fopen(fullfile(outDir, 'Table2_Synthetic_Validation.csv'), 'w');
    fprintf(fid, 'Method,Peak_Optical_Depth,AUC_mean,AUC_std,Pd_at_Pfa_1e-2,Seeds\n');
    for m = 1:nM
        for a = 1:nA
            fprintf(fid, '"%s",%g,%.4f,%.4f,%.4f,%d\n', methods{m, 1}, amps(a), ...
                    T.aucMean(m, a), T.aucStd(m, a), T.pdMean(m, a), nS);
        end
    end
    fclose(fid);
 
    if true
        fig = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1100 440]);
        subplot(1, 2, 1);  hold on;  grid on;  box on;
        for m = 1:nM
            h = errorbar(amps, T.aucMean(m, :), T.aucStd(m, :), '-o');  set(h, 'LineWidth', 1.4);
        end
        set(gca, 'XScale', 'log');  ylim([0.5 1.01]);
        xlabel('Peak injected optical depth (2005 nm)');  ylabel('AUC (mean \pm sd over seeds)');
        legend(methods(:, 1), 'Location', 'southeast');  title('Detection quality vs plume strength');
 
        subplot(1, 2, 2);  hold on;  grid on;  box on;
        D = U.synth('Amplitude', 0.02, 'Seed', seeds(1));
        for m = 1:nM
            [~, ~, info] = methods{m, 2}(D);
            mt = U.metrics(info.raw, D.truth);
            plot(mt.fpr, mt.tpr, 'LineWidth', 1.4);
        end
        plot([0 1], [0 1], 'k:');  xlim([0 0.3]);
        xlabel('False-alarm rate');  ylabel('Probability of detection');
        title('ROC, peak optical depth 0.02');
        U.saveFig(fig, fullfile(outDir, 'Figure5_Synthetic_Validation.png'));
        close(fig);
    end
end
 
