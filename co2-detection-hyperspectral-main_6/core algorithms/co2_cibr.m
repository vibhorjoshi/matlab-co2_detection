function [co2_index, hotspot_mask] = co2_cibr()
    clc; close all;
    mat_file = fullfile(pwd, 'proposed_results.mat');
    load(mat_file, 'cube', 'wavelengths'); wavelengths = wavelengths(:);
    
    [~, b_L] = min(abs(wavelengths - 1980)); 
    [~, b_A] = min(abs(wavelengths - 2005)); 
    [~, b_R] = min(abs(wavelengths - 2030)); 
    
    left = max(cube(:,:,b_L), 0); absorb = max(cube(:,:,b_A), 0); right = max(cube(:,:,b_R), 0);
    co2_index = ((left + right) / 2) - absorb;
    
    co2_index = max(co2_index, 0);
    co2_index = medfilt2(co2_index, [3 3]);
    co2_index = co2_index / (max(co2_index(:)) + 1e-10);
    
    hotspot_mask = bwareaopen(co2_index > prctile(co2_index(:), 95), 10);
    
    figure('Name', 'Stage 1: CIBR', 'Color', 'w');
    subplot(1,2,1); imagesc(co2_index); axis image; colormap('parula'); title('CIBR Map');
    subplot(1,2,2); imshow(hotspot_mask); title('CIBR Hotspots (Scattered Noise)');
end