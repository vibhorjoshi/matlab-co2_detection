function [co2_map, hotspot_mask, info] = co2_sfa(varargin)
%CO2_SFA  Stage 3 - Spectral Fitting Algorithm on continuum-removed spectra.
%   [co2_map, hotspot_mask, info] = co2_sfa(cube, wavelengths)   HxWxB reflectance cube, Bx1 band centres (nm)
%
%   [co2_map, hotspot_mask, info] = co2_sfa(cube, wavelengths)   (same inputs as co2_cibr)

    U = co2_util();
    [D, cfg, ~] = U.prepare(varargin);
    wl = D.wl(:);
    [rows, cols, B] = size(D.cube);  N = rows * cols;
    R2 = reshape(D.cube, N, B);

    [tau, idx] = U.contTau(R2, wl, cfg.features, cfg.nShoulder, D.floor);
    t = U.target(wl(idx), cfg);

    amp = tau * t;                                           % <t, tau>, ||t|| = 1
    switch lower(cfg.sfa.score)
        case 'cosine',    raw = amp ./ (sqrt(sum(tau .^ 2, 2)) + eps);
        case 'amplitude', raw = amp;
        otherwise, error('co2:sfa', 'Unknown cfg.sfa.score "%s".', cfg.sfa.score);
    end

    [co2_map, hotspot_mask, fin] = U.finalize(reshape(raw, rows, cols), D, cfg);
    info = fin;
    info.bands = numel(idx);
end
