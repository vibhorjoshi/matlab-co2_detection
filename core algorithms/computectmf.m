function [mf_image, hotspot_mask, info] = computectmf(varargin)
%COMPUTECTMF  Stage 4 - Cluster-Tuned ACE (CT-ACE), the proposed method.
%   [mf_image, hotspot_mask, info] = computectmf(cube, wavelengths)   HxWxB reflectance cube, Bx1 band centres (nm)
%
%   [mf_image, hotspot_mask, info] = computectmf()
%       No arguments: runs on proposed_results.mat in the repo root. If that
%       file is missing or is not valid reflectance, the built-in synthetic
%       scene (known plume) is used, so the proposed method can ALWAYS be run
%       without the 41 GB flightline.
%   [...] = computectmf('file.mat')            cube + wavelengths
%   [...] = computectmf('f.img', 'f.hdr')      raw ENVI flightline
%   [...] = computectmf(D)                     struct from CO2_LOAD_CUBE
%   [...] = computectmf('synthetic')
%   Options: 'Config', struct | 'Crop', [r1 r2 c1 c2]|'full' | 'Overlay', false
%
%   Pipeline: SWIR band selection (1500-2100 nm, 1800-1950 nm water excluded)
%             -> absorbance -log(R) (natural log, data-driven floor)
%             -> K-means (cfg.ctace.K, fixed seed)
%             -> per-cluster signed ACE (CO2_CT_ACE)
%             -> shared post-processing (CO2_FINALIZE).
%   Setting cfg.ctace.absorbance = false or cfg.ctace.K = 1 gives the
%   ablation variants used in src/run_validation.m.
    U = co2_util();
    [D, cfg, rest] = U.prepare(varargin);
    o = U.opts(struct('Overlay', false), rest);

    sel = U.bands(D.wl, cfg);
    wls = D.wl(sel);  M = numel(wls);
    [rows, cols, ~] = size(D.cube);  N = rows * cols;

    X = reshape(D.cube(:, :, sel), N, M);
    if cfg.ctace.absorbance
        [X, ainfo] = U.absorbance(X, D.floor);
        info.absorbanceFloor = ainfo.floor;
        info.clampedFrac     = ainfo.clampedFrac;
    end
    t = U.target(wls, cfg);
    if ~cfg.ctace.absorbance                       % reflectance domain: absorption is a DROP
        t = -(t - mean(t));  t = t / norm(t);
    end

    [s, cluster, cinfo] = U.ctAce(X, t, wls, cfg);

    [mf_image, hotspot_mask, fin] = U.finalize(reshape(s, rows, cols), D, cfg);
    f = fieldnames(fin);
    for i = 1:numel(f), info.(f{i}) = fin.(f{i}); end
    info.cluster      = reshape(cluster, rows, cols);
    info.clusterSizes = cinfo.sizes;
    info.borrowedGlobal = cinfo.borrowedGlobal;

    if o.Overlay
        U.geoOverlay(hotspot_mask, D.georef, fullfile(U.root(), 'updated_output_file', ...
                        'Figure4_Geospatial_Validation.png'));
    end
end
