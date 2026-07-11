function [sweep, folds] = sgrBifurcation(paramName, values, base)
%SGRBIFURCATION  Equilibrium branches over a parameter sweep.
%   Sweeps one model parameter, solves for every feasible equilibrium at each
%   value with SGREQUILIBRIA, and returns the equilibria grouped by stability.
%   Because each sweep point is solved independently by polynomial root-finding,
%   the branches do not depend on continuation from a previous point and no
%   branch is lost at a fold.
%
%   One function serves every bifurcation diagram in the paper: the recruitment
%   diagram (beta_sg), deterrence (eps), relapse (beta_rg), reintegration
%   (beta_rs) and advocacy (eta) sweeps differ only in PARAMNAME and VALUES.
%
%   Syntax:
%       sweep = SGRBIFURCATION(paramName, values)
%       [sweep, folds] = SGRBIFURCATION(paramName, values, base)
%
%   Description:
%       sweep = SGRBIFURCATION(paramName, values) returns a struct array with
%       one element per swept value, each with fields:
%           value  - the parameter value
%           E      - n-by-3 equilibria [s g r]
%           kinds  - n-by-1 cell array of stability labels
%
%       [sweep, folds] = SGRBIFURCATION(...) also returns the approximate
%       saddle-node (fold) locations: the swept values at which the number of
%       feasible positive equilibria changes.
%
%   Input Arguments:
%       paramName - Name of the parameter to sweep
%           Type: char vector, one of the fields of the parameter struct
%       values    - Values to sweep over
%           Type: numeric vector
%       base      - Baseline parameters (optional; default SGRBASELINE)
%
%   Output Arguments:
%       sweep - 1-by-numel(values) struct array (fields: value, E, kinds)
%       folds - vector of swept values where the positive-equilibrium count changes
%
%   Example:
%       vals = linspace(0.01, 1.05, 400);
%       [sweep, folds] = sgrBifurcation('beta_sg', vals);
%       % folds ~ 0.0757 (saddle-node); the transcritical point is mu + gamma.
%
%   See also SGREQUILIBRIA, SGRMAKEFIGURES.

    if nargin < 3 || isempty(base)
        base = sgrBaseline();
    end
    if ~isfield(base, paramName)
        error('sgrBifurcation:unknownParam', 'Unknown parameter "%s".', paramName);
    end

    n = numel(values);
    sweep = struct('value', cell(1, n), 'E', cell(1, n), 'kinds', cell(1, n));
    nPos = zeros(1, n);

    for i = 1:n
        p = base;
        p.(paramName) = values(i);
        [E, kinds] = sgrEquilibria(p);
        sweep(i).value = values(i);
        sweep(i).E = E;
        sweep(i).kinds = kinds;
        nPos(i) = sum(E(:, 2) > 1e-6);        % count of positive-gang equilibria
    end

    % A fold is where the number of positive equilibria changes.
    folds = [];
    for i = 2:n
        if nPos(i) ~= nPos(i - 1)
            folds(end + 1) = 0.5 * (values(i) + values(i - 1)); %#ok<AGROW>
        end
    end
end
