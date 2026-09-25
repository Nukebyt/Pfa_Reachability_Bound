function prm = KCFAR_Params(c2, varargin)
%KCFAR_PARAMS  Method-of-log-cumulants estimator for K-distributed clutter
%   (single-look, L=1 -- the amplitude products in this project are all
%   single-look, matching every other detector's convention here).
%
%   prm = KCFAR_PARAMS(c2, ...)
%
%   K-distributed amplitude arises from the classic texture model: intensity
%   I | tau ~ Exponential(mean tau)  (single-look),  tau ~ Gamma(a, mu/a),
%   amplitude V = sqrt(I). Shape inversion reuses `_common/inv_trigamma.m`
%   (already relied on by G0CFAR_Params.m) rather than a fresh solver here --
%   same vectorised, domain-safe damped-Newton root of psi(1,.), one fewer
%   place for that iteration to disagree with itself across detectors. Via
%   the mixture representation (no Bessel functions needed -- see
%   KCFAR_TLog.m for why this route was chosen):
%
%       Var(log I) = Var(log tau) + E[Var(log I | tau)]
%                  = psi(1,a) + psi(1,1)          (psi(1,1) = pi^2/6, exact)
%
%   and log(V) = 0.5*log(I), so Var(log V) = 0.25*Var(log I):
%
%       c2 = ( psi(1,a) + pi^2/6 ) / 4     =>     psi(1,a) = 4*c2 - pi^2/6
%
%   so the shape parameter `a` (the K-distribution's order/texture
%   parameter) inverts via the INVERSE TRIGAMMA function. This is exactly
%   CFAR K/legacy/kmolc.m's formula (a = invpsi2(4*c2-psi(1,L),1) at L=1),
%   independently re-derived here from the mixture representation as a
%   cross-check -- both routes agree exactly, which is the reason this
%   derivation is trusted rather than transplanted from the legacy script
%   as-is (see also the header of Main_Estimation_and_Detection_K.m, whose
%   own CFAR step used a DIFFERENT, Gamma-texture-only estimator (nkgmolc)
%   on MAP-filtered data -- not reusable here, since this project's shared
%   front end has no MAP-filtering stage; see KCFAR_Floating.m).
%
%   SUPPORT CONDITION
%   -----------------
%   psi(1,a) -> 0 as a -> infinity and psi(1,a) -> infinity as a -> 0+, so a
%   finite positive solution exists only when  4*c2 - pi^2/6 > 0, i.e.
%   c2 > pi^2/24 (~0.4112). Below that, the implied texture variance would
%   have to be negative -- there is no K-distribution with this little
%   log-domain variance, only a smaller Gamma/Rayleigh-like limit. This is
%   the same kind of support condition BurrCFAR_Params.m and
%   GenGammaCFAR_Params.m already carry for their own models.
%
%   NAME-VALUE OPTIONS
%     'AMin' (0.05), 'AMax' (100) : clamps on the shape parameter `a`.
%
%   OUTPUT (struct prm)
%     a, a_raw, valid, FractionAtMin, FractionAtMax
%
%   See also: KCFAR_TLog, KCFAR_Floating, cfar_front_end

p = inputParser;
addParameter(p, 'AMin', 0.05, @(v) isnumeric(v) && isscalar(v) && v > 0);
addParameter(p, 'AMax', 100,  @(v) isnumeric(v) && isscalar(v) && v > 0);
parse(p, varargin{:});
AMin = p.Results.AMin;
AMax = p.Results.AMax;
if AMax <= AMin
    error('KCFAR_Params:BadClamp', 'AMax must exceed AMin.');
end

PSI11 = pi^2 / 6;           % psi(1,1), exact closed form (L=1)
target = 4*c2 - PSI11;      % = psi(1,a); must be > 0 for a finite a > 0

a_raw = inv_trigamma(target);           % NaN wherever target <= 0 or non-finite
supported = isfinite(a_raw);
a = min(max(a_raw, AMin), AMax);

prm = struct();
prm.Name            = 'K';
prm.a               = a;
prm.a_raw           = a_raw;
prm.valid           = supported & isfinite(a) & a > 0;
prm.FractionInvalid = 1 - mean(prm.valid(:));
prm.FractionAtMin   = mean(a_raw(supported) <= AMin + 1e-9) ;
prm.FractionAtMax   = mean(a_raw(supported) >= AMax - 1e-9);
prm.AMin = AMin;  prm.AMax = AMax;
end
