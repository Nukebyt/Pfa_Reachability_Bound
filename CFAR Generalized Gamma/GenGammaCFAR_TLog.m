function [delta, info] = GenGammaCFAR_TLog(prm, Pfa)
%GENGAMMACFAR_TLOG  Log-domain threshold offset for generalized-gamma clutter.
%
%   [delta, info] = GENGAMMACFAR_TLOG(prm, Pfa)
%
%   Returns delta = T_log - c1 for the decision  x_cut > c1 + delta.
%
%   ---------------------------------------------------------------------
%   DERIVATION -- and how log(sigma) disappears
%   ---------------------------------------------------------------------
%   The original threshold (CFAR Generalized Gamma/legacy, both variants) is
%       T = sigma * ( (1/k) * G )^(1/v),
%       G = gammaincinv(1-Pfa, k)   if v > 0
%           gammaincinv(  Pfa, k)   if v < 0
%   with sigma = exp(c1 - (psi(k) - log k)/v) from the MoLC estimator.
%   Taking logs:
%       log T = log(sigma) + (log G - log k)/v
%             = c1 - (psi(k) - log k)/v + (log G - log k)/v
%             = c1 + (log G - psi(k)) / v
%   The two log(k) terms cancel exactly, leaving
%
%       delta = [ log(gammaincinv(tail, k)) - psi(k) ] / v
%
%   This is the same c1-plus-a-shape-only-offset form every detector in this
%   project reduces to, and it is what makes delta precomputable into a ROM:
%   it depends on (k, v, Pfa) and on nothing about the image's brightness.
%
%   ---------------------------------------------------------------------
%   THE TAIL SWITCH
%   ---------------------------------------------------------------------
%   v < 0 inverts the mapping from the underlying gamma variate to the
%   amplitude, so the upper amplitude tail corresponds to the LOWER gamma
%   tail. Hence gammaincinv(Pfa, k) rather than gammaincinv(1-Pfa, k) there.
%   The original scripts branch on `if powe>0` for exactly this reason.
%   Note the branch is on the CLAMPED v's sign, which is the same as the
%   raw v's sign (clamping only bounds |v|).
%
%   ---------------------------------------------------------------------
%   HARDWARE NOTE (Phase 3/4)
%   ---------------------------------------------------------------------
%   delta = numerator(k, Pfa) / v. Both k and v derive from the same
%   quantised (r, c2) pair, so the whole expression collapses to a single
%   ROM addressed by {Pfa_sel, r_addr, c2_addr} -- or, better, to
%   numerator(k,Pfa) from a 1-D ROM on r multiplied by a reciprocal-|v|
%   obtained from a second 1-D ROM, since v = sign(-c3)*sqrt(psi(1,k)/c2)
%   factorises cleanly into a k-only part and a c2-only part:
%       1/v = sign(-c3) * sqrt(c2) / sqrt(psi(1,k))
%   This is the generalized-gamma analogue of the combined K/C LUT that
%   removed the runtime divider from the Weibull datapath, and it means no
%   divider is needed here either.
%
%   INPUTS
%     prm : struct from GenGammaCFAR_Params (fields k, v, valid)
%     Pfa : scalar, 0 < Pfa < 1
%
%   OUTPUTS
%     delta : HxW map of T_log - c1 (NaN where prm.valid is false)
%     info  : diagnostics -- fraction using each tail branch

if ~isscalar(Pfa) || Pfa <= 0 || Pfa >= 1
    error('GenGammaCFAR_TLog:BadPfa', 'Pfa must be a scalar in (0,1) -- got %g.', Pfa);
end

k = prm.k;
v = prm.v;

% gammaincinv is evaluated only on valid, finite entries; elsewhere delta is
% left NaN and cfar_compose suppresses detection there.
delta = nan(size(k));
ok = prm.valid & isfinite(k) & isfinite(v) & v ~= 0;
if ~any(ok(:))
    info = struct('Pfa', Pfa, 'FractionUpperTail', NaN, 'FractionLowerTail', NaN);
    return;
end

kv = k(ok);
vv = v(ok);

pos = vv > 0;
G = zeros(size(kv));
% v > 0: upper amplitude tail <-> upper gamma tail.
% 'upper' is used rather than gammaincinv(1-Pfa, k) so that small Pfa keeps
% full precision (1-1e-6 would round away most of the tail's information).
if any(pos)
    G(pos) = gammaincinv(Pfa, kv(pos), 'upper');
end
% v < 0: upper amplitude tail <-> LOWER gamma tail.
if any(~pos)
    G(~pos) = gammaincinv(Pfa, kv(~pos), 'lower');
end

d = (log(G) - psi(kv)) ./ vv;
delta(ok) = d;

info = struct();
info.Pfa               = Pfa;
info.FractionUpperTail = mean(pos);
info.FractionLowerTail = mean(~pos);
end
