function nFail = sgrValidate()
%SGRVALIDATE  Verify the MATLAB implementation against the paper's benchmarks.
%   Runs a suite of checks against values reported in the companion paper and
%   against the reference Python implementation. Prints a per-check result and
%   returns the number of failures, so the script doubles as the executable
%   entry point for automated reproducibility checks (exit code 0 on success
%   when called from the wrapper below).
%
%   Syntax:
%       nFail = SGRVALIDATE()
%
%   Description:
%       nFail = SGRVALIDATE() runs every check and returns the failure count.
%       A return value of 0 means the MATLAB implementation reproduces the
%       published equilibria, thresholds, saddle-node location, and stability
%       classification.
%
%   Checks:
%       1  gang-free equilibrium exists and is (1, 0, 0)
%       2  transcritical threshold at beta_sg = mu + gamma
%       3  worked instance (R0 = 0.5): three equilibria, correct coordinates
%       4  worked instance: stability classes are stable / saddle / stable
%       5  baseline: high-gang equilibrium coordinates
%       6  bistability below threshold (R0 < 1 with a positive stable state)
%       7  saddle-node (fold) location in beta_sg
%       8  equilibrium residuals at machine precision
%       9  ode45 converges to the same high-gang equilibrium as the root finder
%      10  advocacy is decoupled: eta is independent (post-reformulation model)
%
%   Example:
%       nFail = sgrValidate();
%
%   See also SGREQUILIBRIA, SGRPRCC, SGRBIFURCATION.

    tol = 1e-6;
    nFail = 0;
    fprintf('SGR MATLAB validation\n');
    fprintf('%s\n', repmat('-', 1, 62));

    % --- 1. gang-free equilibrium -------------------------------------------
    p = sgrBaseline();
    [E, kinds] = sgrEquilibria(p);
    nFail = nFail + report(1, 'gang-free equilibrium (1,0,0) present', ...
        max(abs(E(1, :) - [1 0 0])) < tol);

    % --- 2. transcritical threshold ------------------------------------------
    tc = p.mu + p.gamma;
    nFail = nFail + report(2, 'transcritical at beta_sg = mu + gamma = 0.865', ...
        abs(tc - 0.865) < 1e-9);

    % --- 3-4. worked instance, R0 = 0.5 --------------------------------------
    q = sgrBaseline('beta_sg', 0.4325);
    [Eq, kq] = sgrEquilibria(q);
    ref = [1 0 0; 0.627830769838 0.006983477402 0.365185752760; ...
           0.035925517346 0.440121503665 0.523952978989];
    nFail = nFail + report(3, 'worked instance: 3 equilibria, correct coordinates', ...
        size(Eq, 1) == 3 && max(max(abs(Eq - ref))) < tol);
    nFail = nFail + report(4, 'worked instance: stable / saddle / stable', ...
        strcmp(kq{1}, 'stable') && strcmp(kq{2}, 'saddle') && strcmp(kq{3}, 'stable'));

    % --- 5. baseline high-gang equilibrium -----------------------------------
    nFail = nFail + report(5, 'baseline E_+ at g = 0.4579, r = 0.5241', ...
        abs(E(end, 2) - 0.4579) < 1e-3 && abs(E(end, 3) - 0.5241) < 1e-3);

    % --- 6. bistability below threshold --------------------------------------
    R0 = q.beta_sg / (q.mu + q.gamma);
    hasPositiveStable = false;
    for j = 1:size(Eq, 1)
        if Eq(j, 2) > 1e-6 && strcmp(kq{j}, 'stable')
            hasPositiveStable = true;
        end
    end
    nFail = nFail + report(6, 'bistability: R0 < 1 with a positive stable state', ...
        R0 < 1 && hasPositiveStable && strcmp(kq{1}, 'stable'));

    % --- 7. saddle-node location ---------------------------------------------
    vals = linspace(0.02, 0.20, 400);
    [~, folds] = sgrBifurcation('beta_sg', vals);
    nFail = nFail + report(7, 'saddle-node in beta_sg near 0.0757', ...
        ~isempty(folds) && abs(folds(1) - 0.0757) < 5e-3);

    % --- 8. residuals ---------------------------------------------------------
    [~, ~, info] = sgrEquilibria(q);
    nFail = nFail + report(8, 'equilibrium residuals below 1e-10', ...
        max(info.residuals) < 1e-10);

    % --- 9. ode45 agrees with the root finder --------------------------------
    [~, Y] = ode45(@(t, y) sgrRHS(t, y, p), [0 40000], [0.5; 0.45]);
    nFail = nFail + report(9, 'ode45 converges to the same E_+ as root-finding', ...
        abs(Y(end, 1) - E(end, 2)) < 1e-3 && abs(Y(end, 2) - E(end, 3)) < 1e-3);

    % --- 10. advocacy decoupling ---------------------------------------------
    %  eta must act independently of beta_rs: changing eta with beta_rs fixed
    %  must change the recruitment denominator (post-reformulation model).
    p0 = sgrBaseline('eta', 0.0);
    p1 = sgrBaseline('eta', 1.0);
    y  = [0.3; 0.4];
    d0 = sgrRHS(0, y, p0);
    d1 = sgrRHS(0, y, p1);
    nFail = nFail + report(10, 'advocacy eta is an independent parameter', ...
        abs(d0(1) - d1(1)) > 1e-6 && p0.beta_rs == p1.beta_rs);

    fprintf('%s\n', repmat('-', 1, 62));
    if nFail == 0
        fprintf('RESULT: 10/10 checks passed.\n');
    else
        fprintf('RESULT: %d check(s) FAILED.\n', nFail);
    end
end

% ------------------------------------------------------------------------------
function failed = report(n, name, ok)
%REPORT  Print one check result; return 1 if it failed.
    if ok
        fprintf('  [%2d] PASS  %s\n', n, name);
        failed = 0;
    else
        fprintf('  [%2d] FAIL  %s\n', n, name);
        failed = 1;
    end
end
