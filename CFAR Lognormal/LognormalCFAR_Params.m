function prm = LognormalCFAR_Params(c2, varargin)
%LOGNORMALCFAR_PARAMS  Method-of-log-cumulants estimator for lognormal clutter.
%
%   prm = LOGNORMALCFAR_PARAMS(c2)
%
%   Lognormal clutter model: the AMPLITUDE A is lognormal, i.e.
%       x = log(A) ~ Normal(mu, sigma^2)
%   so the log-cumulants are as trivial as they get:
%       mu       = c1
%       sigma^2  = c2
%   which is exactly what the original CFAR Lognormal/legacy/lgnmolc.m does
%   (it is a two-line pass-through). mu is not handled here at all, because
%   it is carried by the shared additive c1 term -- see LognormalCFAR_TLog.
%
%   This is the ONLY one of the five detectors whose estimator is closed
%   form with no special function, no iteration and no validity condition:
%   every c2 > 0 gives a valid model. That has a direct hardware
%   consequence -- see LognormalCFAR_TLog.m -- and makes Lognormal the
%   natural control case in the Phase 2 comparison.
%
%   INPUTS
%     c2 : HxW map of the log-amplitude variance, from cfar_front_end
%
%   NAME-VALUE OPTIONS
%     'SigmaMin' : lower clamp on sigma (default 0.05)
%     'SigmaMax' : upper clamp on sigma (default 2.0)
%       Clamps bound the threshold a single pathological window can produce.
%       They are applied to sigma (not sigma^2) so the clamp is linear in
%       the quantity that actually multiplies z(Pfa).
%
%   OUTPUT (struct prm)
%     sigma      : HxW CLAMPED standard deviation of log-amplitude
%     sigma_raw  : HxW unclamped sigma (diagnostic, for LUT range sizing)
%     valid      : HxW logical, all true (kept for interface parity with the
%                  other four detectors, which do have real failure modes)
%     FractionAtMin, FractionAtMax, FractionInvalid : scalars
%
%   See also: cfar_front_end, LognormalCFAR_TLog, LognormalCFAR_Floating

p = inputParser;
addParameter(p, 'SigmaMin', 0.05, @(v) isnumeric(v) && isscalar(v) && v > 0);
addParameter(p, 'SigmaMax', 2.0,  @(v) isnumeric(v) && isscalar(v) && v > 0);
parse(p, varargin{:});
sMin = p.Results.SigmaMin;
sMax = p.Results.SigmaMax;
if sMax <= sMin
    error('LognormalCFAR_Params:BadClamp', 'SigmaMax must exceed SigmaMin.');
end

sigma_raw = sqrt(max(c2, 0));
sigma     = min(max(sigma_raw, sMin), sMax);

prm = struct();
prm.Name            = 'Lognormal';
prm.sigma           = sigma;
prm.sigma_raw       = sigma_raw;
prm.valid           = true(size(c2));
prm.FractionAtMin   = mean(sigma_raw(:) <= sMin + 1e-12);
prm.FractionAtMax   = mean(sigma_raw(:) >= sMax - 1e-12);
prm.FractionInvalid = 0;
prm.SigmaMin        = sMin;
prm.SigmaMax        = sMax;
end
