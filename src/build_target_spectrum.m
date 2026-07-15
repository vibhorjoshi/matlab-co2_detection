function sig = build_target_spectrum(wavelengths)
% BUILD_TARGET_SPECTRUM  Dual-band Gaussian CO2 absorption template
%
% Constructs a synthetic CO2 spectral signature for use as the target
% vector d in the matched filter / ACE stages.  The template models the
% two principal CO2 SWIR absorption features:
%
%   Band 1 (near-SWIR)  : amplitude=0.30, centre=1575 nm, sigma=15 nm
%                         (3v1+v3 vibrational combination mode)
%   Band 2 (mid-SWIR)   : amplitude=0.70, centre=2005 nm, sigma=12 nm
%                         (4v3+2delta overtone mode, dominant feature)
%
% The amplitude ratio A2/A1 = 7/3 reflects the relative absorption
% strengths as observed by AVIRIS sensors (Green 2001; Dennison 2013).
%
% The output is L2-normalised so that the matched-filter score is
% expressed in units of signal-to-noise rather than raw energy.
%
% USAGE
%   d = build_target_spectrum(wavelengths)
%
% INPUT
%   wavelengths  [B x 1] or [1 x B]  Band centre wavelengths in nm
%
% OUTPUT
%   sig          [B x 1]  L2-normalised target signature
%
% This is the SINGLE canonical definition used across the project.
% co2_sfa.m and co2_ctmf.m each contain a local copy named
% buildtargetspectrum() so they can run standalone; those copies
% are identical to this function.
%
% REFERENCES
%   Green (2001); Dennison et al. (2013); HITRAN database
% =========================================================================

    wavelengths = wavelengths(:);   % force column vector
    B   = numel(wavelengths);
    sig = zeros(B, 1);

    amp1 = 0.30;  cen1 = 1575;  sig1 = 15;
    amp2 = 0.70;  cen2 = 2005;  sig2 = 12;

    for b = 1:B
        lam    = wavelengths(b);
        sig(b) = amp1 * exp(-0.5*((lam - cen1)/sig1)^2) ...
               + amp2 * exp(-0.5*((lam - cen2)/sig2)^2);
    end

    sig = sig / (norm(sig) + 1e-10);   % L2-normalise
end
