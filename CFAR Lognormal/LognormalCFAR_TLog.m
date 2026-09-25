function [delta, info] = LognormalCFAR_TLog(prm, Pfa)
%LOGNORMALCFAR_TLOG  Log-domain threshold offset for lognormal clutter.
%
%   [delta, info] = LOGNORMALCFAR_TLOG(prm, Pfa)
%
%   Returns delta = T_log - c1, so the detector's decision is
%       x_cut > c1 + delta
%   where x = log(sqrt(I+0.5)). This is the per-distribution "back end" that
%   Phase 3 turns into a ROM; everything before it is shared.
%
%   DERIVATION
%     x = log(A) ~ N(mu, sigma^2),  mu = c1,  sigma = sqrt(c2)
%     P(A > T) = Pfa  <=>  (log T - mu)/sigma = z_{1-Pfa}
%     =>  T_log = log T = c1 + z * sigma
%     =>  delta = z(Pfa) * sigma,     z(Pfa) = sqrt(2)*erfcinv(2*Pfa)
%
%   The erfcinv form is used rather than the algebraically identical
%   norminv(1-Pfa): at Pfa = 1e-6 the argument 2*Pfa = 2e-6 is exact, while
%   1-Pfa rounds to 0.999999 and loses the very precision the tail
%   threshold depends on. (norminv_local.m documents the same point.)
%
%   WHY LOGNORMAL IS THE CHEAPEST DETECTOR IN HARDWARE
%     delta = z * sigma with z a compile-time constant per Pfa. So the whole
%     back end is:  one LUT (c2 -> sqrt(c2)) and one multiply by a constant.
%     There is NO reciprocal and NO divider anywhere -- unlike Weibull
%     (K/C), Generalized Gamma (.../v) and Burr (.../rho), all of which need
%     either a divider or the combined-LUT trick to remove one. Even the
%     single LUT is optional: c2 -> sqrt is cheap enough to compute directly.
%     If Phase 2 shows Lognormal's Pd within a few points of the others, it
%     is very likely the correct production choice on a resource-bound
%     device regardless of what the SAR literature prefers.
%
%   INPUTS
%     prm : struct from LognormalCFAR_Params
%     Pfa : scalar desired probability of false alarm, 0 < Pfa < 0.5
%
%   OUTPUTS
%     delta : HxW map of T_log - c1
%     info  : struct with the scalar constant z used (Phase 3 LUT entry)

if ~isscalar(Pfa) || Pfa <= 0 || Pfa >= 0.5
    error('LognormalCFAR_TLog:BadPfa', ...
        'Pfa must be a scalar in (0, 0.5) -- got %g. (z(Pfa) changes sign at 0.5.)', Pfa);
end

z = sqrt(2) * erfcinv(2 * Pfa);
delta = z .* prm.sigma;

info = struct();
info.z   = z;
info.Pfa = Pfa;
end
