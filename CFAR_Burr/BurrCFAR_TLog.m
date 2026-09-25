function [delta, info] = BurrCFAR_TLog(prm, Pfa)
%BURRCFAR_TLOG  Log-domain threshold offset for Burr-XII clutter.
%
%   [delta, info] = BURRCFAR_TLOG(prm, Pfa)
%
%   Returns delta = T_log - c1 for the decision  x_cut > c1 + delta.
%
%   ---------------------------------------------------------------------
%   DERIVATION
%   ---------------------------------------------------------------------
%   The Burr-XII survival function is S(t) = (1 + (t/eta)^rho)^(-kappa), so
%   S(T) = Pfa gives the legacy threshold (CFAR_Burr/legacy/Main_CFAR_BurrXII.m)
%       T = eta * ( Pfa^(-1/kappa) - 1 )^(1/rho)
%   and the MoLC scale is eta = exp( c1 - (psi(1) - psi(kappa))/rho ).
%   Taking logs, the eta term and the threshold term share the same 1/rho:
%
%       log T = c1 - (psi(1) - psi(kappa))/rho + log(Pfa^(-1/kappa) - 1)/rho
%
%       delta = [ psi(kappa) - psi(1) + log( Pfa^(-1/kappa) - 1 ) ] / rho
%
%   psi(1) = -EulerGamma = -0.5772156649. Once again delta depends only on
%   the shape parameters and Pfa -- never on the clutter level, which is
%   carried entirely by the shared additive c1.
%
%   ---------------------------------------------------------------------
%   THE OVERFLOW THAT MUST BE AVOIDED
%   ---------------------------------------------------------------------
%   Pfa^(-1/kappa) = exp(w) with w = -log(Pfa)/kappa > 0. For Pfa = 1e-6 and
%   a small kappa -- and small kappa is exactly the heavy-tailed regime Burr
%   is chosen FOR -- w runs into the hundreds and exp(w) overflows to Inf
%   long before the model stops being meaningful. kappa = 0.1 at Pfa = 1e-6
%   already gives w = 138, and exp(138) is ~1e60; by kappa = 0.01 it is Inf.
%   Computing log(exp(w) - 1) naively therefore returns Inf, silently killing
%   detection on precisely the windows the detector exists to handle.
%
%   The stable identity used instead is
%       log(e^w - 1) = w + log1p(-e^(-w))
%   which is exact for every w > 0: e^(-w) is in (0,1) so log1p never
%   overflows, and for large w it degrades gracefully to w itself. For small
%   w (large kappa) log1p(-e^-w) carries the correction accurately, so no
%   branch is needed anywhere in the range.
%
%   ---------------------------------------------------------------------
%   HARDWARE NOTE (Phase 3/4)
%   ---------------------------------------------------------------------
%   delta = numerator(kappa, Pfa) / rho, with kappa from a 1-D ROM on the
%   skewness s = c3/c2^1.5 and rho = sqrt((psi(1,kappa)+psi(1,1))/c2). As for
%   the generalized gamma, 1/rho factorises into a kappa-only part and a
%   c2-only part:
%       1/rho = sqrt(c2) / sqrt(psi(1,kappa) + psi(1,1))
%   so the divider disappears into two 1-D ROMs and one multiply -- the same
%   trick as the Weibull combined K/C LUT. Burr remains the most expensive of
%   the five (two shape parameters, two ROMs, a reciprocal-sqrt) which is why
%   it is scheduled last in Phase 5.
%
%   INPUTS
%     prm : struct from BurrCFAR_Params (fields kappa, rho, valid)
%     Pfa : scalar, 0 < Pfa < 1
%
%   OUTPUTS
%     delta : HxW map of T_log - c1 (NaN where prm.valid is false)
%     info  : diagnostics, including the largest w encountered -- worth
%             watching, since it bounds the dynamic range the Phase 3
%             fixed-point format has to cover.

if ~isscalar(Pfa) || Pfa <= 0 || Pfa >= 1
    error('BurrCFAR_TLog:BadPfa', 'Pfa must be a scalar in (0,1) -- got %g.', Pfa);
end

delta = nan(size(prm.kappa));
ok = prm.valid & isfinite(prm.kappa) & isfinite(prm.rho) & prm.kappa > 0 & prm.rho > 0;
if ~any(ok(:))
    info = struct('Pfa', Pfa, 'MaxW', NaN, 'MinW', NaN);
    return;
end

kappa = prm.kappa(ok);
rho   = prm.rho(ok);

w = -log(Pfa) ./ kappa;                  % > 0
log_term = w + log1p(-exp(-w));          % = log(Pfa^(-1/kappa) - 1), overflow-free

delta(ok) = (psi(kappa) - psi(1) + log_term) ./ rho;

info = struct();
info.Pfa  = Pfa;
info.MinW = min(w);
info.MaxW = max(w);
end
