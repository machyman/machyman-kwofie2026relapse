function [E, kinds, info] = sgrEquilibria(p, tol)
%SGREQUILIBRIA  All feasible equilibria of the SGR model, by polynomial elimination.
%   Finds every feasible steady state of the susceptible-gang-reformed model in
%   the homogeneous-contact reduction and classifies each one from the
%   eigenvalues of the Jacobian of the reduced (g, r) system.
%
%   The method is deterministic: the steady-state equations are reduced to a
%   single polynomial in the reformed fraction r, whose roots are found with
%   ROOTS. Unlike a Newton iteration, this recovers every branch, including the
%   saddle equilibrium that separates the basins of attraction, and it does not
%   depend on an initial guess. This matters because the model can be bistable,
%   so a solver seeded from one initial condition can silently miss a branch.
%
%   Derivation:
%       With s = 1 - g - r, the equation dr/dt = 0 solves for g as a rational
%       function of r,
%           g = N(r) / Dn(r),   N(r) = (beta_rs + mu)*r - beta_rs*r^2,
%                               Dn(r) = gamma + (beta_rs - beta_rg)*r,
%       and s = S(r)/Dn(r) with S(r) = beta_rg*r^2 - (beta_rg + k)*r + gamma,
%       where k = gamma + mu. Substituting into dg/dt = 0 (for g > 0) and
%       clearing Dn gives the polynomial
%           P(r) = beta_sg*S(r) - (Dn(r) + eta*r*S(r) + eps*N(r))*(k - beta_rg*r)
%       which is quartic in r when eta > 0 and cubic when eta = 0. Feasible
%       roots satisfy 0 <= r <= 1, g > 0 and s >= 0.
%
%   Syntax:
%       E = SGREQUILIBRIA(p)
%       [E, kinds] = SGREQUILIBRIA(p)
%       [E, kinds, info] = SGREQUILIBRIA(p, tol)
%
%   Description:
%       E = SGREQUILIBRIA(p) returns an n-by-3 array whose rows are the feasible
%       equilibria [s g r], ordered by increasing gang fraction g. The gang-free
%       equilibrium [1 0 0] is always the first row.
%
%       [E, kinds] = SGREQUILIBRIA(p) also returns an n-by-1 cell array of
%       stability labels, each 'stable', 'saddle' or 'unstable'.
%
%       [E, kinds, info] = SGREQUILIBRIA(p, tol) also returns a struct with the
%       polynomial coefficients, the residual of each equilibrium, and the
%       eigenvalues used for classification.
%
%   Input Arguments:
%       p - Model parameters (see SGRBASELINE)
%           Type: struct with fields beta_sg, beta_rg, beta_rs, mu, gamma, eps, eta
%
%       tol - Feasibility and residual tolerance (optional)
%           Type: scalar double
%           Default: 1e-9
%
%   Output Arguments:
%       E     - n-by-3 double, rows [s g r]
%       kinds - n-by-1 cell array of char vectors: 'stable' | 'saddle' | 'unstable'
%       info  - struct with fields: coeffs, residuals, eigenvalues
%
%   Example:
%       p = sgrBaseline('beta_sg', 0.4325);      % worked instance, R0_simp = 0.5
%       [E, kinds] = sgrEquilibria(p);
%       % E_0 = (1,0,0) stable; E_S = (0.628,0.007,0.365) saddle;
%       % E_+ = (0.036,0.440,0.524) stable  -> bistability below threshold
%
%   Note:
%       Uses only base MATLAB (ROOTS, EIG). No toolbox is required, and the
%       function runs unchanged under GNU Octave.
%
%   Reference:
%       Kwofie, Rodriguez Rodriguez, Wang, Hyman and Kang, "Relapse-Driven
%       Bistability in a Nonlinear Model of Gang Dynamics with Deterrence and
%       Reformed-Led Advocacy", MBE-8616.
%
%   See also SGRBASELINE, SGRJACOBIAN, SGRRHS, SGRR0.

    if nargin < 2 || isempty(tol)
        tol = 1e-9;
    end

    a  = p.beta_sg;
    b  = p.beta_rg;
    c  = p.beta_rs;
    m  = p.mu;
    gm = p.gamma;
    ep = p.eps;
    et = p.eta;
    k  = gm + m;

    % Coefficients of the inner factor  Dn + eta*r*S + eps*N  (ascending powers)
    i0 = gm;
    i1 = (c - b) + et * gm + ep * (c + m);
    i2 = -et * (b + k) - ep * c;
    i3 = et * b;

    % P(r) = a*S(r) - (inner)*(k - b*r), expanded (ascending powers of r)
    p0 = a * gm - i0 * k;
    p1 = -a * (b + k) - (i1 * k - i0 * b);
    p2 = a * b - (i2 * k - i1 * b);
    p3 = -(i3 * k - i2 * b);
    p4 = i3 * b;

    coeffs = [p4, p3, p2, p1, p0];          % descending powers, for ROOTS
    while numel(coeffs) > 1 && abs(coeffs(1)) < 1e-14
        coeffs = coeffs(2:end);             % drop leading zeros (eta = 0 -> cubic)
    end

    % The gang-free equilibrium always exists.
    E = [1, 0, 0];

    if numel(coeffs) > 1
        rts = roots(coeffs);
        for j = 1:numel(rts)
            rj = rts(j);
            if abs(imag(rj)) > 1e-8
                continue                     % complex root: not an equilibrium
            end
            r = real(rj);
            if r < -tol || r > 1 + tol
                continue
            end
            Dn = gm + (c - b) * r;
            if abs(Dn) < 1e-12
                continue                     % degenerate: g not determined by r
            end
            N = (c + m) * r - c * r^2;
            g = N / Dn;
            s = 1 - g - r;
            if g < tol || s < -tol
                continue                     % infeasible (need g > 0, s >= 0)
            end
            if isempty(findRow(E, [s, g, r]))
                E = [E; s, g, r];            %#ok<AGROW>
            end
        end
    end

    % Order by gang fraction, then classify.
    [~, ord] = sort(E(:, 2));
    E = E(ord, :);

    n = size(E, 1);
    kinds = cell(n, 1);
    residuals = zeros(n, 1);
    evs = zeros(n, 2);
    for j = 1:n
        g = E(j, 2);
        r = E(j, 3);
        J = sgrJacobian([g; r], p);
        ev = eig(J);
        evs(j, :) = ev(:).';
        re = real(ev);
        if all(re < -1e-9)
            kinds{j} = 'stable';
        elseif any(re > 1e-9) && any(re < -1e-9)
            kinds{j} = 'saddle';
        else
            kinds{j} = 'unstable';
        end
        residuals(j) = norm(sgrRHS(0, [g; r], p), inf);
    end

    info = struct('coeffs', coeffs, 'residuals', residuals, 'eigenvalues', evs);
end

% ------------------------------------------------------------------------------
function idx = findRow(E, row)
%FINDROW  Index of a row of E matching ROW to 1e-6, or empty if none.
    idx = [];
    for j = 1:size(E, 1)
        if max(abs(E(j, :) - row)) < 1e-6
            idx = j;
            return
        end
    end
end
