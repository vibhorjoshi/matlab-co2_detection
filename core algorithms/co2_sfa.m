function [co2_map, binary_map] = co2_sfa(cube, wavelengths)
%CO2_SFA  Spectral Fitting Algorithm – correlation with dual-band CO2 template
%
%   [co2_map, binary_map] = co2_sfa(cube, wavelengths)
%
% Template is a dual Gaussian centred at 1575 nm and 2005 nm (L2-normalised).
% Matching is performed after mean-centring both the observed spectrum and
% the template so that the continuum shape does not dominate the score.

    if nargin < 2
        error('co2_sfa:NotEnoughInputs', 'Provide cube and wavelengths.');
    end

    wavelengths = wavelengths(:);
    [H, W, ~] = size(cube);

    band_mask = wavelengths >= 1500 & wavelengths <= 2100;
    if sum(band_mask) < 10
        band_mask = true(size(wavelengths));
    end
    wl = wavelengths(band_mask);
    X  = double(cube(:,:,band_mask));
    nB = numel(wl);

    % Dual-band CO2 absorption template (identical definition used by CT-ACE)
    amp1 = 0.30; cen1 = 1575; sig1 = 15;
    amp2 = 0.70; cen2 = 2005; sig2 = 12;
    sig = amp1 * exp(-0.5*((wl-cen1)/sig1).^2) + ...
          amp2 * exp(-0.5*((wl-cen2)/sig2).^2);
    sig  = sig / norm(sig);
    ref0 = sig - mean(sig);                 % mean-centred

    % Vectorised per-pixel correlation
    X2d = reshape(X, H*W, nB);
    mu  = mean(X2d, 2);
    sd  = std(X2d, 0, 2);
    sd(sd < 1e-8) = 1;                      % avoid division by zero
    Xn  = (X2d - mu) ./ sd;
    scores = (Xn * ref0) / nB;
    scores = max(scores, 0);

    co2_map = reshape(scores, H, W);
    mx = max(co2_map(:));
    if mx > 0
        co2_map = co2_map / mx;
    end

    thr = prctile(co2_map(:), 95);
    binary_map = bwareaopen(co2_map > thr, 10);
end
