function p = sgrBaseline(varargin)
%SGRBASELINE  Illustrative baseline parameters for the SGR gang-dynamics model.
%   Returns the illustrative (not calibrated) parameter set used for the figures
%   and the sensitivity analysis in Kwofie, Rodriguez Rodriguez, Wang, Hyman and
%   Kang (Mathematical Biosciences and Engineering, MBE-8616), Table 2.
%
%   The model uses the post-decoupling advocacy formulation: advocacy
%   effectiveness ETA is an independent parameter (baseline 0), and does not
%   inherit its scale from the reintegration rate BETA_RS.
%
%   Syntax:
%       p = SGRBASELINE()
%       p = SGRBASELINE('beta_sg', 0.4325)
%
%   Description:
%       p = SGRBASELINE() returns a struct of the illustrative baseline values.
%
%       p = SGRBASELINE(NAME, VALUE, ...) returns the baseline with the named
%       fields overridden, which is the usual way to sweep a parameter.
%
%   Output Arguments:
%       p - Model parameters
%           Type: struct with scalar double fields
%               beta_sg - recruitment rate (week^-1)
%               beta_rg - relapse rate (week^-1)
%               beta_rs - reintegration rate, reformed to susceptible (week^-1)
%               mu      - demographic turnover rate (week^-1)
%               gamma   - reformation rate (week^-1)
%               eps     - perceived-risk deterrence coefficient (dimensionless)
%               eta     - reformed-led advocacy effectiveness (dimensionless)
%
%   Example:
%       p = sgrBaseline();                    % baseline, R0_simp ~ 0.99
%       q = sgrBaseline('beta_sg', 0.4325);   % worked instance, R0_simp = 0.5
%
%   Note:
%       These values are illustrative and chosen to expose the nonlinear
%       bifurcation structure. They are not calibrated to a specific community.
%
%   See also SGREQUILIBRIA, SGRRHS, SGRR0.

    p = struct( ...
        'beta_sg', 0.855, ...   % recruitment (illustrative; R0_simp ~ 0.99)
        'beta_rg', 1.63, ...    % relapse
        'beta_rs', 0.0001, ...  % reintegration (reformed -> susceptible)
        'mu',      0.005, ...   % demographic turnover
        'gamma',   0.86, ...    % reformation
        'eps',     0.95, ...    % perceived-risk deterrence
        'eta',     0.0);        % advocacy effectiveness (independent; baseline off)

    % Name-value overrides (base MATLAB; no inputParser dependency needed)
    if mod(numel(varargin), 2) ~= 0
        error('sgrBaseline:pairs', 'Name-value arguments must come in pairs.');
    end
    for k = 1:2:numel(varargin)
        name = varargin{k};
        if ~isfield(p, name)
            error('sgrBaseline:unknownField', 'Unknown parameter "%s".', name);
        end
        p.(name) = varargin{k + 1};
    end
end
