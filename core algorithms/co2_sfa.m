function [co2_map, binary_map] = co2_sfa()
    clc; close all;
    load(fullfile(pwd, 'proposed_results.mat'), 'cube', 'wavelengths'); 
    [rows, cols, ~] = size(cube); wavelengths = wavelengths(:);
    
    swir_mask = wavelengths >= 1500 & wavelengths <= 2100;
    wl_swir = wavelengths(swir_mask);
    
    % This now calls the local function attached at the bottom of this file
    ref = buildtargetspectrum(wl_swir);
    ref0 = ref - mean(ref);
    
    co2_map = zeros(rows, cols);
    swir_data = double(cube(:,:,swir_mask));
    
    for r = 1:rows
        for c = 1:cols
            spec = squeeze(swir_data(r,c,:));
            if std(spec) > 1e-10
                spec_norm = (spec - mean(spec)) / std(spec);
                co2_map(r,c) = dot(spec_norm, ref0) / numel(wl_swir);
            end
        end
    end
    
    co2_map = max(co2_map, 0);
    co2_map = co2_map / (max(co2_map(:)) + 1e-10);
    binary_map = bwareaopen(co2_map > prctile(co2_map(:), 95), 10);
    
    figure('Name', 'Stage 3: SFA', 'Color', 'w');
    subplot(1,2,1); imagesc(co2_map); axis image; colormap('parula'); title('SFA Map');
    subplot(1,2,2); imshow(binary_map); title('SFA Hotspots');
end

% =========================================================================
% LOCAL FUNCTION: buildtargetspectrum
% =========================================================================
function sig = buildtargetspectrum(wavelengths)
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