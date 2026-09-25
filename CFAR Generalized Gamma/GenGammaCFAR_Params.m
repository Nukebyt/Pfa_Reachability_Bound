function prm = GenGammaCFAR_Params(c2, c3, varargin)
%GENGAMMACFAR_PARAMS  Method-of-log-cumulants estimator for generalized-gamma
%   clutter (shape k, power v, scale sigma).
%
%   prm = GENGAMMACFAR_PARAMS(c2, c3, ...)
%
%   ---------------------------------------------------------------------
%   THE MoLC SYSTEM
%   ---------------------------------------------------------------------
%   For the generalized gamma distribution the log-cumulants are
%       c1 = log(sigma) + (psi(k) - log k)/v
%       c2 = psi(1,k) / v^2
%       c3 = psi(2,k) / v^3
%   (this is what CFAR Generalized Gamma/legacy/gengammamolc.m inverts).
%   Eliminating v between the last two gives a single equation in k alone:
%
%       c3^2 / c2^3  =  psi(2,k)^2 / psi(1,k)^3   =:  g(k)
%
%   ---------------------------------------------------------------------
%   WHY k DEPENDS ONLY ON A SINGLE DIMENSIONLESS RATIO -- and why that
%   matters for the FPGA
%   ---------------------------------------------------------------------
%   r = c3^2/c2^3 = skew^2 is scale-free: it is unchanged by any rescaling
%   of the image. So the shape parameter k is a function of ONE scalar, not
%   of the (c2,c3) pair. Generalized gamma therefore needs a 1-D shape ROM
%   in Phase 3, not a 2-D one -- which is the difference between a few
%   kilobits and a few megabits of on-chip memory. The same argument applies
%   to Burr XII (see BurrCFAR_Params.m); it is the single most important
%   structural fact for making either 3-parameter detector fit a Cyclone V.
%
%   ---------------------------------------------------------------------
%   TWO SOLVERS: 'cubic' (the original) and 'exact'
%   ---------------------------------------------------------------------
%   The original gengammamolc.m does NOT solve g(k) = r. It solves a cubic,
%   via Cardano:
%       a0 k^3 + a1 k^2 + a2 k + a3 = 0
%       a0 = 8c3^2, a1 = 4(3c3^2-2c2^3), a2 = 2(3c3^2-8c2^3), a3 = c3^2-8c2^3
%   Dividing through by c2^3 and collecting in r, that cubic is IDENTICALLY
%       r * (2k+1)^3 = 8 * (k+1)^2      i.e.   r = 8(k+1)^2/(2k+1)^3
%   (verified symbolically and numerically: the cubic's unique positive real
%   root satisfies this exactly). That expression is the LARGE-k ASYMPTOTIC
%   approximation to g(k), obtained from psi(1,k) ~ 1/k + 1/(2k^2) and
%   psi(2,k) ~ -1/k^2 - 1/k^3.
%
%   The approximation is NOT benign at the small k typical of heterogeneous
%   SAR clutter (measured, this file's own relation vs. the exact one):
%       k    = 0.2    0.5    1.0    2.0    5.0    10     50
%       err  = +20%   -4.5%  -8.7%  -5.4%  -1.5%  -0.4%  -0.02%
%   and the two have DIFFERENT domains: g(k) -> 4 as k -> 0+, while the
%   cubic's r(k) -> 8. So for r in (4, 8) the cubic returns a root that the
%   exact MoLC system has no solution for at all.
%
%   Both are provided. 'cubic' is the default so this model reproduces the
%   supplied formula; 'exact' inverts g(k) itself. The Phase 2 comparison
%   runs both and reports the difference, and the choice costs nothing in
%   hardware -- either way Phase 3 stores a 1-D table indexed by r, and the
%   table may as well hold the exact values.
%
%   Neither solver uses fsolve (no Optimization Toolbox installed), and
%   neither uses Cardano's closed form directly: on the three-real-root
%   branch the discriminant goes negative and MATLAB's principal-branch
%   cube root of a complex number does not return the real root, which is a
%   latent bug in the original. Both solvers instead invert the relevant
%   monotone function via mono_table_invert.
%
%   ---------------------------------------------------------------------
%   REMAINING PARAMETERS
%   ---------------------------------------------------------------------
%       v     = sign(-c3) * sqrt(psi(1,k)/c2)
%       sigma = exp(c1 - (psi(k) - log k)/v)
%   v's sign comes from c3's: c3 = psi(2,k)/v^3 and psi(2,k) < 0 always, so
%   v < 0 exactly when c3 > 0. sigma is NOT computed here -- it depends on
%   c1, and c1 is carried by the shared additive term (see
%   GenGammaCFAR_TLog, where the log(sigma) contribution cancels down to
%   -psi(k)/v).
%
%   NAME-VALUE OPTIONS
%     'Solver'  : 'cubic' (default, the original relation) | 'exact'
%     'KMin'    : lower clamp on k (default 0.05)
%     'KMax'    : upper clamp on k (default 50)
%     'VAbsMin' : floor on |v| (default 0.05) -- v appears only as a
%                 divisor in the threshold, so a near-zero v is the one
%                 numerically dangerous case.
%     'VAbsMax' : ceiling on |v| (default 20)
%
%   OUTPUT (struct prm)
%     k, v            : HxW clamped shape and power parameters
%     k_raw, v_raw    : unclamped (diagnostic / LUT range sizing)
%     r               : HxW c3^2/c2^3, the LUT address variable
%     valid           : HxW logical -- false where r is outside the solver's
%                       invertible range (no MoLC solution exists)
%     FractionInvalid, FractionAtMin, FractionAtMax
%
%   See also: cfar_front_end, GenGammaCFAR_TLog, mono_table_invert

p = inputParser;
addParameter(p, 'Solver',  'cubic', @(v) any(strcmpi(v, {'cubic','exact'})));
addParameter(p, 'KMin',    0.05, @(v) isnumeric(v) && isscalar(v) && v > 0);
addParameter(p, 'KMax',    50,   @(v) isnumeric(v) && isscalar(v) && v > 0);
addParameter(p, 'VAbsMin', 0.05, @(v) isnumeric(v) && isscalar(v) && v > 0);
addParameter(p, 'VAbsMax', 20,   @(v) isnumeric(v) && isscalar(v) && v > 0);
parse(p, varargin{:});
solver  = lower(p.Results.Solver);
kMin    = p.Results.KMin;
kMax    = p.Results.KMax;
vAbsMin = p.Results.VAbsMin;
vAbsMax = p.Results.VAbsMax;

%% ---- The scale-free address variable ------------------------------------
r = (c3.^2) ./ (c2.^3);

%% ---- Invert for k -------------------------------------------------------
% The k-grid spans well past the clamps on both sides so that clamping and
% out-of-support are distinguishable: falling off the grid means "no
% solution" (valid=false), hitting a clamp means "solution exists but is
% extreme" (valid=true, saturated).
persistent kgrid
if isempty(kgrid)
    kgrid = logspace(-4, 4, 40001)';
end

switch solver
    case 'cubic'
        % r(k) = 8(k+1)^2/(2k+1)^3 -- strictly decreasing, range (0, 8)
        fh = @(k) 8 * (k + 1).^2 ./ (2*k + 1).^3;
    case 'exact'
        % g(k) = psi(2,k)^2/psi(1,k)^3 -- strictly decreasing, range (0, 4)
        fh = @(k) psi(2, k).^2 ./ psi(1, k).^3;
end

[k_raw, valid] = mono_table_invert(fh, kgrid, r);

%% ---- v from k and c2 ----------------------------------------------------
k = min(max(k_raw, kMin), kMax);

% psi(1,k) > 0 for all k > 0 and c2 > 0, so the sqrt is always real.
v_mag_raw = sqrt(psi(1, k) ./ c2);
v_sign    = sign(-c3);
% c3 == 0 exactly means a perfectly symmetric log-histogram: the generalized
% gamma degenerates (v -> inf, it becomes lognormal). Treat it as invalid
% rather than picking an arbitrary sign.
degenerate = (v_sign == 0) | ~isfinite(v_mag_raw);
v_sign(v_sign == 0) = 1;

v_raw = v_sign .* v_mag_raw;
v     = v_sign .* min(max(v_mag_raw, vAbsMin), vAbsMax);

valid = valid & ~degenerate & isfinite(k_raw) & isfinite(v);

%% ---- Output -------------------------------------------------------------
prm = struct();
prm.Name            = 'GeneralizedGamma';
prm.Solver          = solver;
prm.k               = k;
prm.k_raw           = k_raw;
prm.v               = v;
prm.v_raw           = v_raw;
prm.r               = r;
prm.valid           = valid;
prm.FractionInvalid = 1 - mean(valid(:));
prm.FractionAtMin   = mean(k_raw(valid) <= kMin + 1e-12);
prm.FractionAtMax   = mean(k_raw(valid) >= kMax - 1e-12);
prm.KMin = kMin;  prm.KMax = kMax;
prm.VAbsMin = vAbsMin;  prm.VAbsMax = vAbsMax;
end
