function [co2_index, hotspot_mask] = co2_cibr(cube, wavelengths)
%CO2_CIBR  Continuum Interpolated Band Ratio for CO2 detection (~2005 nm)
%
%   [co2_index, hotspot_mask] = co2_cibr(cube, wavelengths)
%
% Inputs
%   cube         HxWxB hyperspectral cube (reflectance or radiance)
%   wavelengths  Bx1 band centres in nanometres
%
% Outputs
%   co2_index    HxW normalised score map in [0,1]
%   hotspot_mask HxW logical mask (95th percentile + bwareaopen)

    if nargin < 2
        error('co2_cibr:NotEnoughInputs', 'Provide cube and wavelengths.');
    end

    wavelengths = wavelengths(:);
    [H, W, ~] = size(cube);

    % Continuum bands around the strong CO2 feature
    left_wl  = 1980;
    abs_wl   = 2005;
    right_wl = 2030;

    [~, iL] = min(abs(wavelengths - left_wl));
    [~, iA] = min(abs(wavelengths - abs_wl));
    [~, iR] = min(abs(wavelengths - right_wl));

    L = double(cube(:,:,iL));
    A = double(cube(:,:,iA));
    R = double(cube(:,:,iR));

    continuum = (L + R) / 2;
    cibr = continuum - A;
    cibr = max(cibr, 0);                 % rectify negative values

    cibr = medfilt2(cibr, [3 3]);

    mx = max(cibr(:));
    if mx > 0
        co2_index = cibr / mx;
    else
        co2_index = zeros(H, W);
    end

    thr = prctile(co2_index(:), 95);
    hotspot_mask = bwareaopen(co2_index > thr, 10);
end
