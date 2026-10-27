function [gas_map, binary_hotspots] = co2_jrge(cube, wavelengths)
%CO2_JRGE  Joint Reflectance and Gas Estimation (iterative continuum residual)
%
%   [gas_map, binary_hotspots] = co2_jrge(cube, wavelengths)
%
% Operates in the 2000-2100 nm window. Accumulates gas-column estimate
% over three iterations of continuum residual fitting.

    if nargin < 2
        error('co2_jrge:NotEnoughInputs', 'Provide cube and wavelengths.');
    end

    wavelengths = wavelengths(:);
    [H, W, ~] = size(cube);

    band_mask = wavelengths >= 2000 & wavelengths <= 2100;
    if sum(band_mask) < 5
        band_mask = wavelengths >= 1950 & wavelengths <= 2150;
    end
    X = double(cube(:,:,band_mask));
    nBands = size(X, 3);

    gasCol = zeros(H, W);
    nIter  = 3;
    scale  = 0.05;

    for t = 1:nIter
        % Linear continuum between first and last band of the window
        C = zeros(H, W, nBands);
        for b = 1:nBands
            alpha = (b-1) / max(nBands-1, 1);
            C(:,:,b) = (1-alpha)*X(:,:,1) + alpha*X(:,:,end);
        end

        residual = max(C - X, 0);                 % absorption only
        sumRes   = sum(residual, 3);
        tau      = scale * sumRes / nBands;
        gasCol   = gasCol + tau;

        % Remove estimated gas contribution
        shape = residual ./ max(sumRes, 1e-8);
        for b = 1:nBands
            X(:,:,b) = X(:,:,b) + tau .* shape(:,:,b);
        end
    end

    mx = max(gasCol(:));
    if mx > 0
        gas_map = gasCol / mx;
    else
        gas_map = zeros(H, W);
    end

    thr = prctile(gas_map(:), 95);
    binary_hotspots = bwareaopen(gas_map > thr, 10);
end
