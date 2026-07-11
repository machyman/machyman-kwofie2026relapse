function J = sgrJacobian(y, p)
%SGRJACOBIAN  Jacobian of the reduced SGR system at a state.
%   Returns the 2-by-2 Jacobian of the reduced (g, r) system defined in SGRRHS,
%   evaluated at y = [g; r] with s = 1 - g - r. Stability of an equilibrium is
%   determined from the eigenvalues of this matrix: an equilibrium is locally
%   stable when both eigenvalues have negative real part, and a saddle when the
%   real parts have opposite signs.
%
%   The derivatives are computed analytically. Writing D = 1 + eta*r*s + eps*g
%   and noting that s = 1 - g - r depends on both state variables,
%       dD/dg = eta*r*(-1) + eps,     dD/dr = eta*(s - r),
%   which are carried through the quotient rule on the recruitment term.
%
%   Syntax:
%       J = SGRJACOBIAN(y, p)
%
%   Input Arguments:
%       y - State [g; r]
%           Type: 2-by-1 double
%       p - Model parameters (see SGRBASELINE)
%           Type: struct
%
%   Output Arguments:
%       J - Jacobian [d(dg)/dg, d(dg)/dr; d(dr)/dg, d(dr)/dr]
%           Type: 2-by-2 double
%
%   Example:
%       p = sgrBaseline();
%       J = sgrJacobian([0.4579; 0.5241], p);   % at the high-gang equilibrium
%       all(real(eig(J)) < 0)                   % -> true (locally stable)
%
%   See also SGRRHS, SGREQUILIBRIA.

    g = y(1);
    r = y(2);
    s = 1 - g - r;

    D   = 1 + p.eta * r * s + p.eps * g;
    dDg = -p.eta * r + p.eps;          % d/dg of (eta*r*(1-g-r) + eps*g)
    dDr = p.eta * (s - r);             % d/dr of (eta*r*(1-g-r))

    % dg/dt = beta_sg*s*g/D + beta_rg*r*g - (gamma+mu)*g
    %   note ds/dg = -1 and ds/dr = -1
    dgg = p.beta_sg * ((-g + s) * D - s * g * dDg) / D^2 + p.beta_rg * r ...
          - (p.gamma + p.mu);
    dgr = p.beta_sg * ((-g) * D - s * g * dDr) / D^2 + p.beta_rg * g;

    % dr/dt = gamma*g - beta_rg*r*g - beta_rs*r*s - mu*r
    drg = p.gamma - p.beta_rg * r + p.beta_rs * r;      % ds/dg = -1
    drr = -p.beta_rg * g - p.beta_rs * (s - r) - p.mu;  % ds/dr = -1

    J = [dgg, dgr; drg, drr];
end
