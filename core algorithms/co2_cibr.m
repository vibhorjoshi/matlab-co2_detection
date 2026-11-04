function [co2_index, hotspot_mask, info] = co2_cibr(varargin)
%CO2_CIBR  Stage 1 - Continuum Interpolated Band Ratio (three-band).
%   [co2_index, hotspot_mask, info] = co2_cibr(cube, wavelengths)   HxWxB reflectance cube, Bx1 band centres (nm)
%
%   [co2_index, hotspot_mask, info] = co2_cibr()                 fallback data
%   [...] = co2_cibr('file.mat')  |  co2_cibr('f.img','f.hdr')  |  co2_cibr(D)
%   Options: 'Config', struct | 'Crop', [r1 r2 c1 c2] | 'full'
%
%   Band centres (cfg.cibr, nm):  left 1980 | absorption 2005 | right 2030.
%   The continuum at the absorption band is linearly interpolated in
%   wavelength between the shoulders, and the score is the fractional band depth
%       CIBR = 1 - R_A / C_A ,   C_A = R_L + (R_R - R_L) (l_A - l_L)/(l_R - l_L)
%   (positive where the surface absorbs below its local continuum).
%
%   OUTPUT  co2_index    [rows x cols] score in [0,1]
%           hotspot_mask [rows x cols] logical (cfg.threshold rule + bwareaopen)
%           info.raw     signed filtered score (for ROC/AUC), info.bands, ...
    U = co2_util();
    [D, cfg, ~] = U.prepare(varargin);
    wl = D.wl(:);
    [~, iL] = min(abs(wl - cfg.cibr.left));
    [~, iA] = min(abs(wl - cfg.cibr.center));
    [~, iR] = min(abs(wl - cfg.cibr.right));
    if ~(iL < iA && iA < iR)
        error('co2:cibr', 'CIBR needs left < centre < right bands; the cube has %.0f/%.0f/%.0f nm.', ...
              wl(iL), wl(iA), wl(iR));
    end
    wA = (wl(iA) - wl(iL)) / (wl(iR) - wl(iL));
    RL = D.cube(:, :, iL);  RA = D.cube(:, :, iA);  RR = D.cube(:, :, iR);
    CA = RL + (RR - RL) * wA;
    raw = 1 - RA ./ max(CA, D.floor);
 
    [co2_index, hotspot_mask, fin] = U.finalize(raw, D, cfg);
    info = fin;
    info.bands = wl([iL iA iR]);
end
