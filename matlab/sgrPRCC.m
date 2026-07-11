function [prcc, ci, out] = sgrPRCC(N, seed, sampleFile)
%SGRPRCC  PRCC sensitivity of long-term gang prevalence to the model parameters.
%   Computes Partial Rank Correlation Coefficients for the active-gang fraction
%   g* at the high-prevalence equilibrium E_+, using Latin Hypercube Sampling
%   around the illustrative baseline, together with bootstrap confidence
%   intervals. This reproduces Figure 8 and Table 4 of the companion paper.
%
%   The response is evaluated at the equilibrium itself, not by integrating from
%   an initial condition. Because the model is bistable, a response defined by
%   simulation would depend on the initial condition through basin membership;
%   evaluating g* at E_+ makes the sensitivity independent of that choice.
%
%   Sampling ranges follow the paper: five parameters within +/-20% of their
%   illustrative baseline, and advocacy effectiveness eta over [0, 1], since its
%   baseline value is zero and a multiplicative range would be degenerate.
%
%   TOOLBOX INDEPENDENCE. Latin Hypercube sampling, rank transformation with
%   tie averaging, partial correlation, and the bootstrap are all implemented
%   here in base MATLAB. The Statistics and Machine Learning Toolbox functions
%   LHSDESIGN, TIEDRANK, PARTIALCORR, BOOTSTRP and CORR are deliberately NOT
%   used, so this function runs in base MATLAB and under GNU Octave.
%
%   Syntax:
%       prcc = SGRPRCC()
%       [prcc, ci] = SGRPRCC(N, seed)
%       [prcc, ci, out] = SGRPRCC(N, seed, sampleFile)
%
%   Description:
%       prcc = SGRPRCC() runs the default analysis (N = 1000) and returns a
%       6-by-1 vector of PRCC values ordered as
%           [beta_sg; gamma; beta_rg; beta_rs; eps; eta].
%
%       [prcc, ci] = SGRPRCC(N, seed) additionally returns a 6-by-2 array of
%       95% bootstrap percentile confidence intervals.
%
%       [prcc, ci, out] = SGRPRCC(N, seed, sampleFile) reads the Latin Hypercube
%       design from SAMPLEFILE (an N-by-6 CSV) instead of generating one. This is
%       the recommended path for cross-language reproducibility: MATLAB and
%       Python seeded random streams differ, so sharing the design matrix is the
%       only way to obtain identical samples and therefore identical PRCC values.
%
%   Input Arguments:
%       N          - Number of Latin Hypercube samples (default 1000)
%       seed       - Random seed (default 42); ignored when SAMPLEFILE is given
%       sampleFile - Path to an N-by-6 CSV design matrix (optional)
%
%   Output Arguments:
%       prcc - 6-by-1 PRCC values, in PARAM_ORDER
%       ci   - 6-by-2 bootstrap 95% confidence intervals [lo, hi]
%       out  - struct: X (design), y (responses), names, fracPositive, meanG
%
%   Example:
%       % Reproduce the published values exactly by sharing Python's design:
%       [prcc, ci] = sgrPRCC(1000, 42, '../data/prcc_lhs_samples.csv');
%
%   See also SGREQUILIBRIA, SGRBASELINE.

    if nargin < 1 || isempty(N),    N = 1000;  end
    if nargin < 2 || isempty(seed), seed = 42; end
    if nargin < 3, sampleFile = ''; end

    base = sgrBaseline('beta_sg', 0.855);
    names = {'beta_sg', 'gamma', 'beta_rg', 'beta_rs', 'eps', 'eta'};
    k = numel(names);

    % --- sampling bounds: +/-20%, except eta which is absolute [0, 1] ---
    frac = 0.20;
    lo = zeros(1, k);
    hi = zeros(1, k);
    for j = 1:k
        if strcmp(names{j}, 'eta')
            lo(j) = 0.0;
            hi(j) = 1.0;
        else
            v = base.(names{j});
            lo(j) = (1 - frac) * v;
            hi(j) = (1 + frac) * v;
        end
    end

    % --- Latin Hypercube design ---
    if ~isempty(sampleFile) && exist(sampleFile, 'file')
        X = dlmread(sampleFile);            % shared design: exact cross-language parity
        N = size(X, 1);
    else
        X = latinHypercube(N, k, lo, hi, seed);
    end

    % --- response: g* at the high-prevalence equilibrium E_+ ---
    y = zeros(N, 1);
    for i = 1:N
        p = base;
        for j = 1:k
            p.(names{j}) = X(i, j);
        end
        E = sgrEquilibria(p);
        [~, kinds] = sgrEquilibria(p);
        gstar = 0;
        for e = 1:size(E, 1)
            if E(e, 2) > 1e-6 && strcmp(kinds{e}, 'stable') && E(e, 2) > gstar
                gstar = E(e, 2);
            end
        end
        y(i) = gstar;
    end

    % --- PRCC and bootstrap CIs ---
    prcc = partialRankCorr(X, y);
    ci = bootstrapPRCC(X, y, 1000, seed);

    out = struct('X', X, 'y', y, 'names', {names}, ...
                 'fracPositive', mean(y > 1e-6), 'meanG', mean(y));
end

% ------------------------------------------------------------------------------
function X = latinHypercube(N, k, lo, hi, seed)
%LATINHYPERCUBE  Stratified Latin Hypercube design (replaces LHSDESIGN).
%   Each column is partitioned into N equal-probability strata; one uniform
%   draw is taken per stratum and the strata are then randomly permuted, which
%   is the defining property of a Latin Hypercube design.
    rand('twister', seed);   %#ok<RAND>  % Octave/MATLAB-portable seeding
    X = zeros(N, k);
    for j = 1:k
        strata = ((0:N-1)' + rand(N, 1)) / N;    % one draw inside each stratum
        X(:, j) = lo(j) + (hi(j) - lo(j)) * strata(randperm(N));
    end
end

% ------------------------------------------------------------------------------
function R = rankTransform(x)
%RANKTRANSFORM  Ranks with averaged ties (replaces TIEDRANK).
    n = numel(x);
    [sorted, idx] = sort(x(:));
    R = zeros(n, 1);
    i = 1;
    while i <= n
        j = i;
        while j < n && sorted(j + 1) == sorted(i)
            j = j + 1;                       % extent of the tie block
        end
        R(idx(i:j)) = (i + j) / 2;           % average rank within the tie block
        i = j + 1;
    end
end

% ------------------------------------------------------------------------------
function prcc = partialRankCorr(X, y)
%PARTIALRANKCORR  PRCC of each column of X with y (replaces PARTIALCORR/CORR).
%   Rank-transforms every variable, then for parameter i correlates the residuals
%   of rank(X_i) and rank(y) after regressing each on the remaining ranked
%   parameters (with an intercept).
    [n, k] = size(X);
    Rx = zeros(n, k);
    for j = 1:k
        Rx(:, j) = rankTransform(X(:, j));
    end
    Ry = rankTransform(y);

    prcc = zeros(k, 1);
    for i = 1:k
        others = setdiff(1:k, i);
        Z = [ones(n, 1), Rx(:, others)];
        rx = Rx(:, i) - Z * (Z \ Rx(:, i));   % residuals after removing the others
        ry = Ry - Z * (Z \ Ry);
        sx = std(rx);
        sy = std(ry);
        if sx > 1e-12 && sy > 1e-12
            prcc(i) = mean((rx - mean(rx)) .* (ry - mean(ry))) / (sx * sy) ...
                      * (n / (n - 1));        % match the population/sample scaling
            % guard against round-off outside [-1, 1]
            prcc(i) = max(min(prcc(i), 1), -1);
        else
            prcc(i) = 0;
        end
    end
end

% ------------------------------------------------------------------------------
function ci = bootstrapPRCC(X, y, nBoot, seed)
%BOOTSTRAPPRCC  Percentile bootstrap CIs for the PRCCs (replaces BOOTSTRP).
    rand('twister', seed + 1);   %#ok<RAND>
    [n, k] = size(X);
    boot = zeros(nBoot, k);
    for b = 1:nBoot
        idx = randi(n, n, 1);
        boot(b, :) = partialRankCorr(X(idx, :), y(idx)).';
    end
    ci = zeros(k, 2);
    for j = 1:k
        col = sort(boot(:, j));
        ci(j, 1) = col(max(1, round(0.025 * nBoot)));
        ci(j, 2) = col(min(nBoot, round(0.975 * nBoot)));
    end
end
