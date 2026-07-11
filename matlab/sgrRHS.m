function dy = sgrRHS(~, y, p)
%SGRRHS  Right-hand side of the reduced SGR gang-dynamics model.
%   Evaluates the two-dimensional reduction of the susceptible-gang-reformed
%   model on the simplex s + g + r = 1, so that s is recovered as s = 1 - g - r
%   and only (g, r) are integrated. The recruitment term is reduced by
%   perceived-risk deterrence and by reformed-led advocacy through the
%   denominator D = 1 + eta*r*s + eps*g.
%
%   Governing equations (with s = 1 - g - r):
%       dg/dt = beta_sg*s*g/D + beta_rg*r*g - (gamma + mu)*g
%       dr/dt = gamma*g - beta_rg*r*g - beta_rs*r*s - mu*r
%
%   Syntax:
%       dy = SGRRHS(t, y, p)
%
%   Description:
%       dy = SGRRHS(t, y, p) returns the time derivative [dg; dr] at state
%       y = [g; r] for the parameters p. The signature is ODE-solver ready:
%       the time argument is ignored (the system is autonomous), so the
%       function can be passed directly to ODE45.
%
%   Input Arguments:
%       t - Time (ignored; the system is autonomous)
%       y - State [g; r]
%           Type: 2-by-1 double, with g, r >= 0 and g + r <= 1
%       p - Model parameters (see SGRBASELINE)
%           Type: struct
%
%   Output Arguments:
%       dy - Time derivative [dg; dr]
%           Type: 2-by-1 double
%
%   Example:
%       p = sgrBaseline();
%       [t, y] = ode45(@(t, y) sgrRHS(t, y, p), [0 2000], [0.5; 0.45]);
%       plot(t, y(:, 1)); xlabel('time (weeks)'); ylabel('g');
%
%   Note:
%       The bifurcation diagrams and the sensitivity analysis do not integrate
%       the model: equilibria are obtained directly by root-finding (see
%       SGREQUILIBRIA). This function exists for trajectory illustrations and
%       for independent cross-checks of the equilibria.
%
%   See also SGREQUILIBRIA, SGRJACOBIAN, SGRBASELINE.

    g = y(1);
    r = y(2);
    s = 1 - g - r;

    D = 1 + p.eta * r * s + p.eps * g;

    dg = p.beta_sg * s * g / D + p.beta_rg * r * g - (p.gamma + p.mu) * g;
    dr = p.gamma * g - p.beta_rg * r * g - p.beta_rs * r * s - p.mu * r;

    dy = [dg; dr];
end
