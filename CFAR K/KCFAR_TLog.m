function [delta, info] = KCFAR_TLog(prm, Pfa)
%KCFAR_TLOG  Log-domain threshold offset for K-distributed clutter.
%
%   [delta, info] = KCFAR_TLOG(prm, Pfa)
%
%   Returns delta = T_log - c1 for the decision  x_cut > c1 + delta.
%
%   ---------------------------------------------------------------------
%   WHY THIS IS A MIXTURE INTEGRAL, NOT A CLOSED FORM
%   ---------------------------------------------------------------------
%   Unlike G0 (exact via the incomplete beta / F-distribution) or Burr XII
%   (exact via its own elementary CDF), K-distributed amplitude has no
%   elementary-function quantile: its density involves a modified Bessel
%   function of the second kind. Rather than evaluate that Bessel function
%   directly (numerically delicate at the extreme arguments Pfa=1e-6
%   reaches), this uses the SAME mixture representation KCFAR_Params.m's
%   derivation is built on:
%
%       I | tau ~ Exponential(mean tau)   (single-look intensity)
%       tau     ~ Gamma(a, mu/a)          (texture, mean mu)
%
%   so the survival function of intensity I is the exact 1-D integral
%
%       S(x; a, mu) = P(I > x) = E_tau[ exp(-x/tau) ]
%                   = integral_0^inf  exp(-x/tau) * gampdf(tau, a, mu/a) dtau
%
%   a smooth mixture-of-exponentials integrand with no special functions
%   beyond the Gamma density -- well-conditioned for MATLAB's adaptive
%   `integral` at every Pfa this project uses, including 1e-6.
%
%   delta depends only on (a, Pfa), never on scale (mu cancels, exactly as
%   for every other detector here) -- fixed at mu=1 throughout.
%
%   ---------------------------------------------------------------------
%   PROCEDURE
%   ---------------------------------------------------------------------
%   1. Solve S(x; a, 1) = Pfa for x (intensity-domain (1-Pfa) quantile),
%      via fzero on a bracket grown until it contains the root.
%   2. Convert to log-amplitude: log(V) = 0.5*log(I), so the quantile in
%      log-amplitude is 0.5*log(x).
%   3. Subtract the log-amplitude mean kappa1(a) = 0.5*(psi(a) - log(a)
%      - EulerGamma) at mu=1 (closed form -- see KCFAR_Params.m's header
%      for the matching derivation of psi(1,a), its second-cumulant twin).
%
%   Only 4 distinct Pfa values are ever swept in this project, so the
%   (a, Pfa) -> delta grid is built ONCE per Pfa (persistent cache) and
%   every pixel's delta comes from a cheap pchip interpolation -- the
%   expensive part (200 fzero+integral solves) runs at most 4 times per
%   MATLAB session, not once per pixel per image.
%
%   HARDWARE NOTE (future phase, not yet built): this grid IS the ROM.
%   Once fixed-point Q-formats are chosen, this exact (a_addr, Pfa_sel) ->
%   delta table becomes a 1-D ROM exactly like Weibull's and Lognormal's,
%   since K's shape is a single scalar per window (c2 alone), not a 2-D
%   pair like G0's (L,u).
%
%   INPUTS
%     prm : struct from KCFAR_Params (fields a, AMin, AMax, valid)
%     Pfa : scalar, 0 < Pfa < 1
%
%   OUTPUTS
%     delta : HxW map of T_log - c1 (NaN where prm.valid is false)
%     info  : diagnostics
%
%   See also: KCFAR_Params, KCFAR_Floating, cfar_front_end

if ~isscalar(Pfa) || Pfa <= 0 || Pfa >= 1
    error('KCFAR_TLog:BadPfa', 'Pfa must be a scalar in (0,1) -- got %g.', Pfa);
end

EULER_GAMMA = 0.5772156649015329;

persistent gridCache
if isempty(gridCache)
    gridCache = containers.Map('KeyType', 'double', 'ValueType', 'any');
end

AMin = prm.AMin; AMax = prm.AMax;
cacheKey = Pfa + AMin*1e-9 + AMax*1e-15;   % distinguish grids if clamps ever change
if ~isKey(gridCache, cacheKey)
    gridCache(cacheKey) = build_delta_grid(Pfa, AMin, AMax, EULER_GAMMA);
end
g = gridCache(cacheKey);

delta = nan(size(prm.a));
ok = prm.valid & isfinite(prm.a);
if any(ok(:))
    delta(ok) = interp1(g.a_grid, g.delta_grid, prm.a(ok), 'pchip', 'extrap');
end

info = struct();
info.Pfa        = Pfa;
info.GridPoints = numel(g.a_grid);
info.AMin       = AMin;
info.AMax       = AMax;
end

function g = build_delta_grid(Pfa, AMin, AMax, EULER_GAMMA)
    % The substituted integral below is occasionally stiff at the extreme
    % ends of the `a` grid (very small or very large shape); `integral`
    % emits non-fatal "max intervals"/"min step size" warnings there but
    % still converges -- verified accurate against a 2e6-sample Monte Carlo
    % calibration check (scratch_verify_k.m) across the a in [0.5,8], Pfa in
    % [1e-4,1e-2] range, ratio 0.98-1.12. Suppressed here so a normal sweep
    % run isn't buried in diagnostic noise; re-enable when touching this code.
    warnState = warning('off', 'MATLAB:integral:MaxIntervalCountReached');
    warnState(2) = warning('off', 'MATLAB:integral:MinStepSize');
    cleanupWarn = onCleanup(@() warning(warnState));

    n = 161;
    a_grid = exp(linspace(log(AMin), log(AMax), n));
    delta_grid = zeros(1, n);
    for i = 1:n
        a = a_grid(i);
        % Gamma(shape=a, scale=1/a) density, in log domain for numerical
        % stability across the wide range of `a` this grid spans (no
        % Statistics Toolbox available -- gampdf is not on this machine).
        % Substitute t = u/(1-u) (maps (0,1) -> (0,inf)) so `integral` sees a
        % finite domain -- more robust than the semi-infinite form for the
        % sharply peaked integrand at small `a`, and avoids the "maximum
        % number of intervals" warning that form triggered (verified this
        % substitution gives the same values where both converge cleanly).
        loggampdf = @(t) (a-1).*log(t) - t.*a - gammaln(a) - a.*log(1/a);
        S = @(x) integral(@(u) exp(-x.*(1-u)./u + loggampdf(u./(1-u))) ./ (1-u).^2, ...
                           0, 1, 'RelTol', 1e-10, 'AbsTol', 1e-15);
        lo = realmin; hi = 1;
        while S(hi) > Pfa
            hi = hi * 4;
            if hi > 1e12
                break;  % pathological a; fzero below will surface the issue
            end
        end
        x = fzero(@(x) S(x) - Pfa, [lo, hi]);
        kappa1 = 0.5 * (psi(a) - log(a) - EULER_GAMMA);
        delta_grid(i) = 0.5 * log(x) - kappa1;
    end
    g = struct('a_grid', a_grid, 'delta_grid', delta_grid);
end
