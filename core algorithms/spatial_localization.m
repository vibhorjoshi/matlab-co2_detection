% =========================================================================
% HYPERSPECTRAL CO2 DETECTION PIPELINE & SPATIAL VALIDATION FRAMEWORK
% =========================================================================
clc; clear; close all;

%% ── CONFIGURATION ────────────────────────────────────────────────────────
BIN_FILE  = 'f250923t01p00r13_rfl.bin';   
HDR_FILE  = 'f250923t01p00r13_rfl.hdr';

% Scene dimensions (from HDR)
N_SAMP  = 1937;      % samples (columns)
N_LINES = 24068;     % lines   (rows)
N_BANDS = 224;
PIXEL_M = 14.4;      % metres per pixel
UL_E    = 577561.59; % UTM Easting of upper-left pixel centre
UL_N    = 4228899.20;% UTM Northing of upper-left pixel centre
UTM_ZONE = 11;       % UTM zone number
INTERLEAVE = 'bil';  % Band Interleaved by Line

% Crop window (rows = lines, cols = samples)
ROW_START = 10000;   ROW_END = 11000;  % 1000 lines
COL_START =   500;   COL_END =  1500;  % 1000 samples
DS = 4;              % downsampling factor → 250×250 working array

% Algorithm parameters
% Tikhonov regularisation
LAMBDA    = 1e-6;    
% CTMF clusters
K_MEANS   = 4;       
% JRGE smoothing spline
SPLINE_P  = 0.90;    
% JRGE iteration step
ALPHA_JR  = 0.05;    
% JRGE iterations
T_ITER    = 3;       
rng(1);

% Band indices (from actual HDR wavelength vector)
B_RED   = 29;   B_GRN = 19;   B_BLU = 8;    % visible RGB
B_SW_R  = 194;  B_SW_G = 137; B_SW_B = 89;  % SWIR false colour
B_CIBR_L = 175; B_CIBR_A = 179; B_CIBR_R = 183; % CIBR bands
CO2_SWIR_START = 128; CO2_SWIR_END = 185;   % 1562–2110 nm for SFA

% P95 threshold for final analysis
P95_THR = 95;

%% ── LOAD AVIRIS DATA (BIL format) ───────────────────────────────────────
fprintf('=== Loading AVIRIS BIL data ===\n');
fprintf('  File    : %s\n', BIN_FILE);
fprintf('  Crop    : rows %d–%d, cols %d–%d\n', ...
    ROW_START,ROW_END,COL_START,COL_END);
assert(isfile(BIN_FILE), 'Binary file not found: %s\nSet BIN_FILE path.', BIN_FILE);

nR_crop = ROW_END  - ROW_START;   % 1000 lines
nC_crop = COL_END  - COL_START;   % 1000 samples

% multibandread for BIL: rows=lines, cols=samples
cube_raw = multibandread(BIN_FILE, [N_LINES N_SAMP N_BANDS], ...
    'float32',  0, INTERLEAVE, 'ieee-le', ...
    {'Row','Range',    [ROW_START ROW_END-1]}, ...
    {'Column','Range', [COL_START COL_END-1]});
cube = double(cube_raw);

% Scale if stored as integers (data type 4 = float32, already scaled)
if max(cube(:)) > 10
    cube = cube / 10000;
end

% Mask invalid
cube(cube <= 0 | cube > 2 | ~isfinite(cube)) = NaN;

% Spatial downsample
cube = cube(1:DS:end, 1:DS:end, :);
[rows, cols, bands] = size(cube);
fprintf('  Working array: %d × %d × %d (DS=%d)\n', rows, cols, bands, DS);

%% ── EXTRACT WAVELENGTH VECTOR FROM HDR ──────────────────────────────────
fprintf('Reading wavelength vector from HDR ...\n');
txt = fileread(HDR_FILE);
tok = regexp(txt, 'wavelength\s*=\s*\{([^}]+)\}', 'tokens');
wl_str = strtrim(strsplit(tok{1}{1}, {',',' ',char(10),char(13)}));
wavelength = str2double(wl_str(~cellfun(@isempty, wl_str)));
wavelength(isnan(wavelength)) = [];
wavelength = wavelength(:);
fprintf('  Bands: %d, range: %.1f – %.1f nm\n', ...
    numel(wavelength), wavelength(1), wavelength(end));

%% ── GEOREFERENCING ───────────────────────────────────────────────────────
% UL corner of the crop window (pixel centre of first retained pixel)
crop_UL_E = UL_E + (COL_START - 1) * PIXEL_M;
crop_UL_N = UL_N - (ROW_START - 1) * PIXEL_M;
pix_m_ds  = PIXEL_M * DS;   % 57.6 m after DS=4

% World limits in UTM metres
xW = [crop_UL_E,  crop_UL_E + cols * pix_m_ds];
yW = [crop_UL_N - rows * pix_m_ds,  crop_UL_N];
fprintf('  Crop UL: E=%.1f m, N=%.1f m\n', crop_UL_E, crop_UL_N);
fprintf('  Pixel spacing after DS: %.1f m\n', pix_m_ds);

% UTM → geographic conversion (requires Mapping Toolbox)
use_mapping_tb = license('test','MAP_Toolbox');
if use_mapping_tb
    utmZ = sprintf('%dN', UTM_ZONE);
    [X_grid, Y_grid] = meshgrid( ...
        linspace(xW(1), xW(2), cols), ...
        linspace(yW(2), yW(1), rows));
    [lat_grid, lon_grid] = utm2ll(X_grid, Y_grid, UTM_ZONE);
    fprintf('  Approx centre: Lat=%.3f, Lon=%.3f\n', ...
        mean(lat_grid(:)), mean(lon_grid(:)));
else
    warning('Mapping Toolbox not available — skipping geo-overlay figures.');
end

%% ── BUILD IMAGE CHANNELS ─────────────────────────────────────────────────
pStr = @(x,p1,p2) max(0,min(1,(x-prctile(x(:),p1))./(prctile(x(:),p2)-prctile(x(:),p1)+1e-8)));

% True-colour RGB
R_ch = pStr(cube(:,:,B_RED), 2, 98);
G_ch = pStr(cube(:,:,B_GRN), 2, 98);
B_ch = pStr(cube(:,:,B_BLU), 2, 98);
rgbImg = cat(3, R_ch, G_ch, B_ch);

% SWIR false-colour (R=2200nm, G=1650nm, B=1200nm)
SW_R = pStr(cube(:,:,B_SW_R), 2, 98);
SW_G = pStr(cube(:,:,B_SW_G), 2, 98);
SW_B = pStr(cube(:,:,B_SW_B), 2, 98);
swirImg = cat(3, SW_R, SW_G, SW_B);

%% ════════════════════════════════════════════════════════════════════════
%% STAGE 1: CIBR
%% ════════════════════════════════════════════════════════════════════════
fprintf('\n[1/4] CIBR ...\n');
R_L = cube(:,:,B_CIBR_L);
R_A = cube(:,:,B_CIBR_A);
R_R = cube(:,:,B_CIBR_R);
R_L(isnan(R_L)) = 0;
R_A(isnan(R_A)) = 0;
R_R(isnan(R_R)) = 0;
continuum  = (R_L + R_R) / 2;
cibr_raw   = continuum - R_A;
cibr_raw(isnan(cibr_raw)) = 0;
p2  = prctile(cibr_raw(:), 2);
p98 = prctile(cibr_raw(:), 98);
cibrScore = (cibr_raw - p2) / (p98 - p2 + 1e-10);
cibrScore  = max(0, min(1, cibrScore));
cibrScore  = medfilt2(cibrScore, [3 3]);
tau_cibr   = graythresh(cibrScore);
cibrMask   = cibrScore > tau_cibr;
fprintf('  Coverage: %.1f%%\n', 100*sum(cibrMask(:))/numel(cibrMask));

%% ════════════════════════════════════════════════════════════════════════
%% STAGE 2: JRGE
%% ════════════════════════════════════════════════════════════════════════
fprintf('[2/4] JRGE (per-pixel spline, %d iterations) ...\n', T_ITER);
absB  = CO2_SWIR_START:CO2_SWIR_END;
wl_idx = 1:bands;
pix2d    = reshape(cube, rows*cols, bands);
g_den    = zeros(rows*cols, 1);
for i = 1:rows*cols
    spec = pix2d(i,:);
    if any(isnan(spec)) || std(spec,'omitnan') < 1e-10
        continue;
    end
    sp  = csaps(wl_idx, spec, SPLINE_P);
    r_s = fnval(sp, wl_idx);
    for t = 1:T_ITER
        r_s = r_s + ALPHA_JR*(spec - r_s);
    end
    g_est = median(r_s(absB)) - median(spec(absB));
    g_den(i) = max(g_est, 0);
    if mod(i, 5000) == 0
        fprintf('  %.0f%%\n', 100*i/(rows*cols));
    end
end
jrgeMap   = reshape(g_den, rows, cols);
g_mn      = min(jrgeMap(:));  g_mx = max(jrgeMap(:));
jrgeScore = (jrgeMap - g_mn) / (g_mx - g_mn + 1e-10);
jrgeScore(isnan(jrgeScore)) = 0;
tau_jrge   = graythresh(jrgeScore);
jrgeMask   = jrgeScore > tau_jrge;
fprintf('  Coverage: %.1f%%\n', 100*sum(jrgeMask(:))/numel(jrgeMask));

%% ════════════════════════════════════════════════════════════════════════
%% STAGE 3: SFA
%% ════════════════════════════════════════════════════════════════════════
fprintf('[3/4] SFA (dual-band 1562–2110 nm) ...\n');
swirIdx  = CO2_SWIR_START:CO2_SWIR_END;
wl_swir  = wavelength(swirIdx);
N_swir   = numel(swirIdx);
A1=0.30; c1=1575; s1=15;
A2=0.70; c2=2005; s2=12;
d_sfa = A1*exp(-0.5*((wl_swir-c1)/s1).^2) + A2*exp(-0.5*((wl_swir-c2)/s2).^2);
d_sfa = d_sfa / norm(d_sfa);
d0    = d_sfa - mean(d_sfa);
sfaMap = zeros(rows, cols);
swirCube = cube(:,:,swirIdx);
for r = 1:rows
    for c = 1:cols
        x = double(squeeze(swirCube(r,c,:)));
        if std(x) < 1e-10, continue; end
        x_s = (x - mean(x)) / std(x);
        sfaMap(r,c) = dot(x_s, d0) / N_swir;
    end
end
sfaMap(sfaMap<0) = 0;
sfaScore = sfaMap / (max(sfaMap(:)) + 1e-10);
tau_sfa   = graythresh(sfaScore);
sfaMask   = sfaScore > tau_sfa;
fprintf('  Coverage: %.1f%%\n', 100*sum(sfaMask(:))/numel(sfaMask));

%% ════════════════════════════════════════════════════════════════════════
%% STAGE 4: CTMF
%% ════════════════════════════════════════════════════════════════════════
fprintf('[4/4] CTMF (k=%d, Tikhonov λ=%.0e) ...\n', K_MEANS, LAMBDA);
d_ctmf = zeros(bands, 1);
for b = 1:bands
    lam = wavelength(b);
    d_ctmf(b) = A1*exp(-0.5*((lam-c1)/s1)^2) + A2*exp(-0.5*((lam-c2)/s2)^2);
end
d_ctmf = d_ctmf / norm(d_ctmf);
sfaW    = 1 + sfaScore;    % SFA-informed weighting
pixFull = reshape(cube, rows*cols, bands);
for b = 1:bands
    pixFull(:,b) = pixFull(:,b) .* sfaW(:);
end
for b = 1:bands
    col_b = pixFull(:,b);
    bad   = isnan(col_b);
    if any(bad), pixFull(bad,b) = nanmedian(col_b); end
end
[clIdx,~] = kmeans(pixFull, K_MEANS, 'MaxIter',200,'Replicates',5,'Distance','sqeuclidean');
mf_out = zeros(rows*cols,1);
for k = 1:K_MEANS
    mem = (clIdx == k);
    if sum(mem) < bands+5, continue; end
    X_k  = pixFull(mem,:);
    mu_k = mean(X_k,1)';
    C_k  = cov(X_k) + LAMBDA*eye(bands);
    w_k  = C_k \ d_ctmf;
    mf_out(mem) = (X_k - mu_k') * w_k;
end
ctmfRaw   = reshape(mf_out, rows, cols);
ctmfRaw(ctmfRaw<0) = 0;
mn_c = min(ctmfRaw(:)); mx_c = max(ctmfRaw(:));
ctmfScore = (ctmfRaw - mn_c) / (mx_c - mn_c + 1e-10);
ctmfScore(isnan(ctmfScore)) = 0;

% Otsu binary mask
tau_otsu  = graythresh(ctmfScore);
ctmfMask  = ctmfScore > tau_otsu;

% P95 binary mask (for connected-component analysis)
tau_p95   = prctile(ctmfScore(:), P95_THR);
maskP95   = ctmfScore > tau_p95;
fprintf('  Coverage (Otsu): %.1f%% | Coverage (P95): %.1f%%\n', ...
    100*sum(ctmfMask(:))/numel(ctmfMask), ...
    100*sum(maskP95(:))/numel(maskP95));

%% ─── CONNECTED COMPONENTS ON P95 MASK ────────────────────────────────────
CC   = bwconncomp(maskP95);
stats = regionprops(CC, 'Area','Centroid','Perimeter','BoundingBox');
areas = [stats.Area];
[~, sortedIdx] = sort(areas,'descend');

% Build colour label map (largest = red, rest = blue)
labelMap = zeros(rows, cols, 'uint8');
for ii = 1:min(CC.NumObjects, 50)
    idx_cc = sortedIdx(ii);
    if ii == 1
        labelMap(CC.PixelIdxList{idx_cc}) = 2;   % red = largest
    else
        labelMap(CC.PixelIdxList{idx_cc}) = 1;   % blue = others
    end
end
LCC_pct = 100 * areas(sortedIdx(1)) / (rows*cols);
compactness = sum(areas) / (CC.NumObjects * pi * (sqrt(areas(sortedIdx(1))/pi))^2);
fprintf('  Connected components (P95): %d | Largest: %.2f%% | Compactness: %.3f\n', ...
    CC.NumObjects, LCC_pct, min(compactness,1));

% Centroid of largest component → geographic position
cen = stats(sortedIdx(1)).Centroid;  % [col, row]
cen_E = crop_UL_E + (cen(1)*DS - 1) * PIXEL_M;
cen_N = crop_UL_N - (cen(2)*DS - 1) * PIXEL_M;
fprintf('  Largest component UTM centre: E=%.1f m, N=%.1f m\n', cen_E, cen_N);

%% ─── MAPPING TOOLBOX FIGURES (Real AVIRIS Image Overlay) ─────────────────
if use_mapping_tb
    fprintf('\nGenerating High-Detail Real Image overlay ...\n');
    
    fig_sat = figure('Name','Detailed Real Map','Units','centimeters', ...
        'Position',[2 2 28 20],'Color','white');
    
    axM = axes(fig_sat);
    
    % 1. Create Latitude and Longitude vectors for the axes
    lon_vec = linspace(min(lon_grid(:)), max(lon_grid(:)), cols);
    lat_vec = linspace(max(lat_grid(:)), min(lat_grid(:)), rows);
    
    % 2. Plot the REAL AVIRIS true-color image as the background
    imagesc(axM, lon_vec, lat_vec, rgbImg);
    set(axM, 'YDir', 'normal'); % Keep north pointing up
    hold(axM, 'on');
    
    % 3. Extract latitude and longitude for the P95 hotspot pixels
    lat_hotspot = lat_grid(maskP95);
    lon_hotspot = lon_grid(maskP95);
    
    % 4. Plot the plume as a dense grid of semi-transparent red dots
    scatter(axM, lon_hotspot, lat_hotspot, 12, 'r', 'filled', 'MarkerFaceAlpha', 0.4);
    
    % 5. Calculate centroids for callout boxes
    cen_lat = mean(lat_hotspot(:));
    cen_lon = mean(lon_hotspot(:));
    
    % Find an extreme western point for a second callout (matching your image)
    [~, min_lon_idx] = min(lon_hotspot);
    west_lat = lat_hotspot(min_lon_idx);
    west_lon = lon_hotspot(min_lon_idx);
    
    % 6. Add Callout Text Boxes matching the uploaded image style
    text(axM, cen_lon, cen_lat, sprintf('Latitude %.4f\nLongitude %.3f', cen_lat, cen_lon), ...
        'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 3, 'FontSize', 9, 'FontWeight', 'bold');
        
    text(axM, west_lon, west_lat, sprintf('Latitude %.4f\nLongitude %.3f', west_lat, west_lon), ...
        'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 3, 'FontSize', 9, 'FontWeight', 'bold');
    
    % Formatting the axes to look like a geographic map
    xlabel(axM, 'Longitude', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel(axM, 'Latitude', 'FontSize', 11, 'FontWeight', 'bold');
    title(axM, '\textbf{High-Resolution Overlay on Real AVIRIS Imagery}', ...
        'Interpreter', 'latex', 'FontSize', 13);
    
    grid(axM, 'on');
    set(axM, 'Layer', 'top'); % Put grid lines over the image
    
    % 7. Export the graphic (This will now be instant and crash-free)
    try
        exportgraphics(fig_sat, 'geo_fig_detailed_real_image.png', 'Resolution', 600);
    catch
        print(fig_sat, 'geo_fig_detailed_real_image.png', '-dpng', '-r600');
    end
    
    fprintf('  Saved detailed real image figure.\n');
end
%% ── LOCAL HELPER FUNCTIONS ───────────────────────────────────────────────
function cmap = viridis_colormap(n)
    % Approximate viridis colormap
    t  = linspace(0,1,n)';
    cmap = [0.267+0.733*t, 0.005+0.874*t.^0.6, 0.329*(1-t)+0.002*t];
    cmap = max(0,min(1,cmap));
end

function [lat,lon] = utm2ll(E,N,zone)
    % Converts UTM Easting/Northing to Latitude/Longitude (WGS-84)
    k0=0.9996; a=6378137; e2=0.00669438;
    x = E - 500000;  y = N;
    M = y/k0;
    mu= M/(a*(1-e2/4-3*e2^2/64));
    n1= (1-sqrt(1-e2))/(1+sqrt(1-e2));
    lat_r = mu+(3/2*n1-27/32*n1^3).*sin(2*mu)+(21/16*n1^2).*sin(4*mu);
    N1= a./sqrt(1-e2.*sin(lat_r).^2);
    T1= tan(lat_r).^2;
    C1= e2./(1-e2).*cos(lat_r).^2;
    R1= a*(1-e2)./(1-e2.*sin(lat_r).^2).^1.5;
    D = x./(N1*k0);
    lat_r = lat_r - N1.*tan(lat_r)./R1.*(D.^2/2 -(5+3*T1+10*C1-4*C1.^2)*D.^4/24);
    lat = rad2deg(lat_r);
    lon_0 = (zone-1)*6-180+3;
    lon = lon_0 + rad2deg(D-( 1+2*T1+C1).*D.^3/6 ./ cos(lat_r));
end