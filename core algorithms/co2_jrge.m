function [gas_map, hotspot_mask, info] = co2_jrge(varargin)
%CO2_JRGE  Stage 2 - Joint Reflectance and Gas Estimation.
%   [gas_map, hotspot_mask, info] = co2_jrge(cube, wavelengths)   HxWxB reflectance cube, Bx1 band centres (nm)
%
%   [gas_map, hotspot_mask, info] = co2_jrge(cube, wavelengths)   (same inputs as co2_cibr)
%
%   Beer-Lambert model in absorbance space, A(l) = -log R(l):
%       A(l) = b(l) + g * t(l)
%   b : smooth surface background (polynomial of degree cfg.jrge.degree)
%   t : CO2 template, g : gas column (template-normalised units, signed)
%   The background is fitted ONLY on bands outside +/- cfg.jrge.guard nm of the
%   template centres (so the gas feature cannot be absorbed by the fit); the
%   gas column is the template-weighted absorbance residual
%       g = <t, A - b> / <t, t> .
%   cfg.jrge.iterations refits b after removing the current gas estimate.
%   Fully vectorised (chunked); needs no Curve Fitting Toolbox.
    U = co2_util();
    [D, cfg, ~] = U.prepare(varargin);
    j   = cfg.jrge;
    sel = U.bands(D.wl, cfg);
    wls = D.wl(sel);  M = numel(wls);
    [rows, cols, ~] = size(D.cube);  N = rows * cols;

    R2 = reshape(D.cube(:, :, sel), N, M);
    A  = U.absorbance(R2, D.floor);
    t  = U.target(wls, cfg);

    x  = (wls - mean(wls)) / (0.5 * (max(wls) - min(wls)));
    B  = x(:) .^ (0:j.degree);                              % M x (deg+1)
    sh = true(M, 1);
    for k = 1:numel(cfg.template.center)
        sh = sh & abs(wls - cfg.template.center(k)) > j.guard;
    end
    if nnz(sh) < j.degree + 3
        error('co2:jrge', 'Too few shoulder bands (%d) for a degree-%d background.', nnz(sh), j.degree);
    end
    P  = pinv(B(sh, :));                                    % (deg+1) x nShoulder

    g = zeros(N, 1);
    for s0 = 1:j.chunk:N
        ii = s0:min(N, s0 + j.chunk - 1);
        Ai = A(ii, :);  gi = zeros(numel(ii), 1);
        for it = 1:j.iterations
            Y    = Ai - gi * t';
            coef = Y(:, sh) * P';                           % N x (deg+1)
            gi   = ((Ai - coef * B') * t) / (t' * t);
        end
        g(ii) = gi;
    end

    [gas_map, hotspot_mask, fin] = U.finalize(reshape(g, rows, cols), D, cfg);
    info = fin;
    info.shoulderBands = nnz(sh);
    info.fractionPositive = mean(fin.raw(:) > 0);
end
