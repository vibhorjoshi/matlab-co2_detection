function [mf_image, binaryMap] = computectmf(varargin)
%COMPUTECTMF  Cluster-Tuned Adaptive Coherence Estimator (CT-ACE)
%
% Recommended calls
%   [mf, mask] = computectmf()                          % uses proposed_results.mat
%   [mf, mask] = computectmf(cube, wavelengths)         % in-memory cube
%   [mf, mask] = computectmf(imgFile, hdrFile)          % raw ENVI pair
%
% The function implements:
%   1. Absorbance transform (-log R) with adaptive floor
%   2. Water-vapour exclusion (1800-1950 nm)
%   3. K-means (K=4) + per-cluster signed ACE
%   4. Rectification, median filter, 95th-percentile threshold, bwareaopen

    % -------------------- Parse inputs --------------------
    cube = [];
    wavelengths = [];

    if nargin == 0
        % No-argument path promised by README
        matFile = locate_fallback();
        S = load(matFile, 'cube', 'wavelengths');
        cube = S.cube;
        wavelengths = S.wavelengths(:);
    elseif nargin == 2 && (ischar(varargin{1}) || isstring(varargin{1}))
        imgFile = char(varargin{1});
        hdrFile = char(varargin{2});
        if exist('hypercube','file') == 2
            hcube = hypercube(imgFile, hdrFile);
            cube  = hcube.DataCube;
            wavelengths = hcube.Wavelength(:);
        else
            info = enviinfo(hdrFile);
            cube = multibandread(imgFile, ...
                [info.Height, info.Width, info.Bands], ...
                'float32', 0, 'bsq', 'ieee-le');
            wavelengths = info.Wavelength(:);
        end
    elseif nargin >= 2 && isnumeric(varargin{1})
        cube = varargin{1};
        wavelengths = varargin{2}(:);
    else
        error('computectmf:InvalidInput', ...
            'Call as computectmf(), computectmf(cube,wl), or computectmf(img,hdr).');
    end

    [H, W, B] = size(cube);
    wavelengths = wavelengths(:);
    if numel(wavelengths) ~= B
        error('computectmf:SizeMismatch', 'Wavelength length must equal number of bands.');
    end

    % -------------------- Scale handling --------------------
    % Some shipped cubes have very small magnitudes (~1e-5).
    % Rescale so that median reflectance lies near 0.3 (typical surface value)
    % before taking the logarithm. This prevents mass clamping.
    medVal = median(cube(:), 'omitnan');
    if medVal > 0 && medVal < 1e-3
        scaleFactor = 0.3 / medVal;
        cube = cube * scaleFactor;
    end

    % -------------------- 1. Band selection + absorbance --------------------
    valid = (wavelengths >= 1500 & wavelengths <= 2100) & ...
            ~(wavelengths >= 1800 & wavelengths <= 1950);
    if sum(valid) < 10
        warning('computectmf:FewBands', 'Few valid bands; using full cube.');
        valid = true(B,1);
    end
    wl_valid = wavelengths(valid);
    R = double(cube(:,:,valid));

    % Adaptive floor: use 1 % of the positive data instead of a hard 1e-5
    pos = R(R > 0);
    if isempty(pos)
        floorVal = 1e-5;
    else
        floorVal = max(prctile(pos, 1), 1e-8);
    end
    A_cube = -log(max(R, floorVal));          % optical depth
    nValid = sum(valid);

    % Dual-band CO2 target signature (identical to SFA)
    amp1=0.30; cen1=1575; sig1=15;
    amp2=0.70; cen2=2005; sig2=12;
    co2_sig = amp1*exp(-0.5*((wl_valid-cen1)/sig1).^2) + ...
              amp2*exp(-0.5*((wl_valid-cen2)/sig2).^2);
    co2_sig = co2_sig / norm(co2_sig);
    co2_sig = co2_sig(:);

    % -------------------- 2. K-means clustering (K=4) --------------------
    reshapedA = reshape(A_cube, H*W, nValid);
    rng(1);                                   % fixed seed for reproducibility
    [clusterIdx, ~] = kmeans(reshapedA, 4, 'MaxIter', 200, ...
                             'Distance', 'sqeuclidean', 'Replicates', 1);

    % -------------------- 3. Per-cluster ACE --------------------
    mf_vec = zeros(H*W, 1);
    for k = 1:4
        idx = find(clusterIdx == k);
        if numel(idx) < nValid + 5
            continue;                         % too few pixels for stable covariance
        end
        Xk = reshapedA(idx, :);
        mu = mean(Xk, 1);
        Xc = Xk - mu;

        Ck = cov(Xk) + eye(nValid)*1e-6;      % regularised
        w  = Ck \ co2_sig;                    % backslash (stable)

        num   = Xc * w;
        Xc_w  = (Ck \ Xc')';
        den1  = sum(Xc .* Xc_w, 2);
        den2  = co2_sig' * w;
        score = sign(num) .* (num.^2) ./ (den1 .* den2 + 1e-8);

        mf_vec(idx) = score;
    end

    mf_image = reshape(mf_vec, H, W);

    % -------------------- 4. Post-processing --------------------
    mf_image = max(mf_image, 0);
    mf_image = medfilt2(mf_image, [3 3]);
    mx = max(mf_image(:));
    if mx > 0
        mf_image = mf_image / mx;
    end

    thr = prctile(mf_image(:), 95);
    binaryMap = bwareaopen(mf_image > thr, 10);
end

% -------------------------------------------------------------------------
function matFile = locate_fallback()
% Locate proposed_results.mat relative to this function or current folder
    candidates = {
        fullfile(fileparts(mfilename('fullpath')), '..', 'proposed_results.mat')
        fullfile(pwd, 'proposed_results.mat')
        fullfile(pwd, '..', 'proposed_results.mat')
        };
    for i = 1:numel(candidates)
        if exist(candidates{i}, 'file')
            matFile = candidates{i};
            return;
        end
    end
    error('computectmf:NoFallback', ...
        ['proposed_results.mat not found. Place it in the repository root ', ...
         'or call computectmf(cube, wavelengths).']);
end      
