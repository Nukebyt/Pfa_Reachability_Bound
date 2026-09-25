function [delta, info] = G0CFAR_TLog(prm, Pfa)
%G0CFAR_TLOG  Log-domain threshold offset for G0 clutter.
%
%   [delta, info] = G0CFAR_TLOG(prm, Pfa)
%
%   Returns delta = T_log - c1 for the decision  x_cut > c1 + delta.
%
%   ---------------------------------------------------------------------
%   DERIVATION -- and a cancellation that removes two logs and the scale
%   ---------------------------------------------------------------------
%   For G0-distributed amplitude, -alpha*X^2/gamma_scale ~ F(2L, -2*alpha),
%   so the legacy threshold (CFAR_G0/legacy/Main_CFAR_G0.m) is
%       T = sqrt( (-gamma_scale/alpha) * finv(1-Pfa, 2L, -2*alpha) )
%   with, from the MoLC estimator,
%       log(gamma_scale) = 2*c1 - psi(L) + psi(u) + log(L),      u = -alpha.
%   Writing Fq for the F quantile and taking logs:
%       log T = 0.5*[ log(gamma_scale) - log(u) + log(Fq) ]
%             = c1 + 0.5*[ -psi(L) + psi(u) + log L - log u + log Fq ]
%
%   Now express Fq through the incomplete beta. With d1 = 2L, d2 = 2u,
%       x  = I^{-1}_{1-Pfa}(L, u)          (regularised incomplete beta inverse)
%       Fq = (d2/d1) * x/(1-x) = (u/L) * x/(1-x)
%   so log Fq = log u - log L + log(x/(1-x)) and BOTH log u and log L cancel:
%
%       delta = 0.5 * [ psi(u) - psi(L) + log( x / (1-x) ) ]
%
%   The scale parameter gamma_scale never has to be formed at all -- which
%   also means the estimator never needs c1, confirming the shared-front-end
%   structure: the entire G0 back end is a function of (L, u, Pfa) only.
%
%   ---------------------------------------------------------------------
%   WHY betaincinv AND NOT finv
%   ---------------------------------------------------------------------
%   finv needs the Statistics Toolbox, which is not installed. betaincinv is
%   a core specfun built-in. The substitution above is exact, not an
%   approximation.
%
%   NUMERICAL POINT: for Pfa = 1e-6, x is within 1e-6 of 1, so computing
%   1-x by subtraction would destroy every significant digit of the very
%   quantity the tail threshold depends on. The beta symmetry
%       I_x(L,u) = 1 - I_{1-x}(u,L)
%   gives 1-x DIRECTLY and accurately as betaincinv(Pfa, u, L) -- note the
%   swapped parameter order and the small first argument. log(x/(1-x)) is
%   then formed as log1p(-y) - log(y) with y = 1-x, which is accurate for
%   every Pfa in range.
%
%   ---------------------------------------------------------------------
%   HARDWARE NOTE (Phase 3/4)
%   ---------------------------------------------------------------------
%   delta depends only on (L, u, Pfa). In 'L1' mode L is a compile-time
%   constant and u comes from a 1-D ROM on c2, so delta collapses to a
%   single ROM addressed by {Pfa_sel, c2_addr} -- no divider, no beta
%   function, nothing transcendental at run time. That is the G0 variant to
%   build first. In 'LA' mode (L, u) is a genuine 2-D shape pair, the only
%   one of the five detectors that does not reduce to one address variable;
%   it needs either a 2-D ROM or an on-chip solver, and its cost should be
%   weighed against the Phase 2 accuracy difference between the two modes.
%
%   INPUTS
%     prm : struct from G0CFAR_Params (fields u, L, valid)
%     Pfa : scalar, 0 < Pfa < 1
%
%   OUTPUTS
%     delta : HxW map of T_log - c1 (NaN where prm.valid is false)
%     info  : diagnostics

if ~isscalar(Pfa) || Pfa <= 0 || Pfa >= 1
    error('G0CFAR_TLog:BadPfa', 'Pfa must be a scalar in (0,1) -- got %g.', Pfa);
end

delta = nan(size(prm.u));
ok = prm.valid & isfinite(prm.u) & isfinite(prm.L) & prm.u > 0 & prm.L > 0;
if ~any(ok(:))
    info = struct('Pfa', Pfa, 'MinBetaTail', NaN, 'MaxBetaTail', NaN);
    return;
end

u = prm.u(ok);
L = prm.L(ok);

% y = 1 - x, computed directly in the accurate direction (see header).
y = betaincinv(Pfa, u, L);
y = min(max(y, realmin), 1 - eps);      % keep both logs finite at the extremes

log_odds = log1p(-y) - log(y);          % = log(x/(1-x))
delta(ok) = 0.5 * (psi(u) - psi(L) + log_odds);

info = struct();
info.Pfa         = Pfa;
info.MinBetaTail = min(y);
info.MaxBetaTail = max(y);
end
