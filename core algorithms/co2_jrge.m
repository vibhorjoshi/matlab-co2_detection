function [gas_map, binary_hotspots] = co2_jrge()
    clc; close all;
    load(fullfile(pwd, 'proposed_results.mat'), 'cube', 'wavelengths'); 
    [rows, cols, ~] = size(cube); wavelengths = wavelengths(:);
    
    co2Mask = wavelengths >= 2000 & wavelengths <= 2100;
    wlW = wavelengths(co2Mask); nBw = sum(co2Mask);
    wlN = (wlW - wlW(1)) / (wlW(end) - wlW(1)); wlN = wlN(:)'; 
    
    X = double(reshape(cube(:,:,co2Mask), rows*cols, nBw));
    gasCol = zeros(rows*cols, 1);
    
    for t = 1:3
        slope = X(:,end) - X(:,1); C = X(:,1) + slope .* wlN; 
        residual = max(C - X, 0); % Fixes the all-zero output bug
        tau = 0.05 * sum(residual, 2) / nBw; gasCol = gasCol + tau;
        shape = residual ./ (sum(residual, 2) + 1e-10); 
        X = X + tau .* shape; 
    end
    
    gas_map = reshape(gasCol, rows, cols);
    gas_map = (gas_map - min(gas_map(:))) / (max(gas_map(:)) - min(gas_map(:)) + 1e-10);
    binary_hotspots = bwareaopen(gas_map > prctile(gas_map(:), 95), 10);
    
    figure('Name', 'Stage 2: JRGE', 'Color', 'w');
    subplot(1,2,1); imagesc(gas_map); axis image; colormap('parula'); title('JRGE Map');
    subplot(1,2,2); imshow(binary_hotspots); title('JRGE Hotspots');
end