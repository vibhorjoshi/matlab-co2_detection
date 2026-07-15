function [mf_image, binaryMap] = co2_ctmf()
    clc; close all;
    load(fullfile(pwd, 'proposed_results.mat'), 'cube', 'wavelengths'); 
    [rows, cols, ~] = size(cube); wavelengths = wavelengths(:);
    
    % Proper Preprocessing: Water Masking & Absorbance Transform
    valid_bands = (wavelengths >= 1500 & wavelengths <= 2100) & ~(wavelengths >= 1800 & wavelengths <= 1950);
    num_valid = sum(valid_bands);
    A_cube = -log(max(double(cube(:,:,valid_bands)), 1e-5));
    
    co2_sig = buildtargetspectrum(wavelengths(valid_bands));
    
    % Clustering
    reshapedA = reshape(A_cube, rows*cols, num_valid);
    rng(42); [clusterIdx, ~] = kmeans(reshapedA, 4, 'MaxIter', 100, 'Distance', 'sqeuclidean');
    
    % CT-ACE Localization (Fixes the "negligible improvement" issue)
    ace_output = zeros(rows*cols, 1);
    for k = 1:4
        idx = (clusterIdx == k);
        X_c = reshapedA(idx, :) - mean(reshapedA(idx, :), 1);
        if sum(idx) < num_valid + 5, continue; end
        
        invC = inv(cov(reshapedA(idx, :)) + eye(num_valid)*1e-5);
        w = invC * co2_sig;
        
        num = X_c * w; 
        den1 = sum(X_c .* (X_c * invC), 2); 
        den2 = co2_sig' * w; 
        ace_output(idx) = sign(num) .* (num.^2) ./ (den1 * den2 + 1e-8);
    end
    
    mf_image = reshape(ace_output, rows, cols);
    mf_image = max(mf_image, 0); 
    mf_image = medfilt2(mf_image, [3 3]);
    mf_image = mf_image / (max(mf_image(:)) + 1e-10);
    
    binaryMap = bwareaopen(mf_image > prctile(mf_image(:), 95), 10);
    
    figure('Name', 'Stage 4: Proposed CTMF (CT-ACE)', 'Color', 'w');
    subplot(1,2,1); imagesc(mf_image); axis image; colormap('parula'); title('Proposed Score Map');
    subplot(1,2,2); imshow(binaryMap); title('Cohesive Plume Isolation');
    
    % Trigger Geospatial Visualization
    mapgeospatial_overlay(binaryMap);
end
