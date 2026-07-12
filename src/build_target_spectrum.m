function sig = build_target_spectrum(wavelengths)
    % Builds a dual-band Gaussian template for CO2 (1575 nm & 2005 nm)
    amp1 = 0.30; cen1 = 1575; sig1 = 15;
    amp2 = 0.70; cen2 = 2005; sig2 = 12;
    
    sig = zeros(length(wavelengths), 1);
    for b = 1:length(wavelengths)
        lam = wavelengths(b);
        sig(b) = amp1*exp(-0.5*((lam-cen1)/sig1)^2) + amp2*exp(-0.5*((lam-cen2)/sig2)^2);
    end
    sig = sig / norm(sig);
end