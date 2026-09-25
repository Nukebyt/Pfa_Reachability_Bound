function [detection_map, threshold_map, prm, stats, logdomain] = ...
    GenGammaCFAR_Floating(I, sli, guard, Pfa, varargin)
%GENGAMMACFAR_FLOATING  Single-pass windowed Generalized-Gamma CFAR detector.
%
%   Floating-point reference model for the streaming FPGA implementation,
%   with the same contract as WeibullCFAR_Floating.m. Composed from
%       cfar_front_end        -- shared with all five detectors
%       GenGammaCFAR_Params   -- (c2,c3) -> k, v  via the scale-free r = c3^2/c2^3
%       GenGammaCFAR_TLog     -- delta = [log gammaincinv(.,k) - psi(k)] / v
%
%   This is one of the two detectors slated for FPGA implementation first.
%
%   DIFFERENCES FROM CFAR Generalized Gamma/legacy/Main_CFAR_Generalized_Gamma*.m
%   ---------------------------------------------------------------------
%   1. Vectorised single pass instead of the per-pixel double loop.
%   2. The legacy CA-CFAR normalisation (`condi = mean + std*thres`) is
%      dropped -- see LognormalCFAR_Floating.m for the full reasoning; it is
%      applied uniformly across all five detectors so Phase 2 compares like
%      with like.
%   3. Cardano's closed form is replaced by a monotone inversion of the same
%      relation. Cardano is not merely slower here, it is wrong on part of
%      the domain: when the cubic's discriminant goes negative (three real
%      roots) the original's `((-q/2)+tmp)^(1/3)` takes MATLAB's principal
%      complex cube root, which does not reconstruct the real root.
%   4. The validity test is kept but restated. The original guards with
%      `3*c3^2 <= 8*c2^3`, i.e. r <= 8/3, and writes a zero (no detection)
%      when it fails; the `_nocond` variant drops the guard and takes abs(k)
%      instead. Here the condition is the solver's actual invertible range
%      (r < 8 for 'cubic', r < 4 for 'exact'), which is the real
%      out-of-support boundary, and windows outside it are counted in
%      stats.FractionInvalid rather than quietly zeroed. The original's
%      stricter r <= 8/3 corresponds to k >~ 0.4 and is available as
%      'RMax', 8/3 if the legacy behaviour is wanted exactly.
%
%   NAME-VALUE OPTIONS
%     'Solver'  : 'cubic' (default, matches the supplied formula) | 'exact'
%     'KMin' (0.05), 'KMax' (50), 'VAbsMin' (0.05), 'VAbsMax' (20)
%     'RMax'    : optional extra upper bound on r = c3^2/c2^3, applied on top
%                 of the solver's own range. Pass 8/3 to reproduce the
%                 legacy `3*c3^2 <= 8*c2^3` guard exactly. Default Inf.
%     'FrontEnd', 'FrontEndOpts', 'Verbose'   -- as for the other detectors
%
%   See also: cfar_front_end, cfar_compose, GenGammaCFAR_Params, GenGammaCFAR_TLog

p = inputParser;
addParameter(p, 'Solver',       'cubic');
addParameter(p, 'KMin',         0.05);
addParameter(p, 'KMax',         50);
addParameter(p, 'VAbsMin',      0.05);
addParameter(p, 'VAbsMax',      20);
addParameter(p, 'RMax',         Inf, @(v) isnumeric(v) && isscalar(v) && v > 0);
addParameter(p, 'FrontEnd',     [],  @(v) isempty(v) || isstruct(v));
addParameter(p, 'FrontEndOpts', {},  @iscell);
addParameter(p, 'Verbose',      false, @(v) islogical(v) && isscalar(v));
parse(p, varargin{:});

if isempty(p.Results.FrontEnd)
    feArg = p.Results.FrontEndOpts;
else
    feArg = p.Results.FrontEnd;
end
rMax = p.Results.RMax;

paramsFcn = @(c2, c3, skew) apply_rmax( ...
    GenGammaCFAR_Params(c2, c3, ...
        'Solver',  p.Results.Solver, ...
        'KMin',    p.Results.KMin,    'KMax',    p.Results.KMax, ...
        'VAbsMin', p.Results.VAbsMin, 'VAbsMax', p.Results.VAbsMax), ...
    rMax);
tlogFcn = @(prm, pf) GenGammaCFAR_TLog(prm, pf);

[detection_map, threshold_map, prm, stats, logdomain] = ...
    cfar_compose(I, sli, guard, Pfa, paramsFcn, tlogFcn, 'GeneralizedGamma', feArg);

if p.Results.Verbose
    v = prm.valid;
    fprintf('GenGammaCFAR_Floating: %dx%d, sli=%d guard=%d N=%d Pfa=%g solver=%s\n', ...
        size(I,1), size(I,2), sli, guard, stats.NumReferenceCells, Pfa, prm.Solver);
    fprintf('  k : min=%.4f mean=%.4f max=%.4f | v : mean=%.4f (%.1f%% negative)\n', ...
        min(prm.k(v)), mean(prm.k(v)), max(prm.k(v)), mean(prm.v(v)), 100*mean(prm.v(v) < 0));
    fprintf('  r = c3^2/c2^3 : median=%.4f  p99=%.4f\n', ...
        median(prm.r(v)), prctile(prm.r(v), 99));
    fprintf('  INVALID windows (no MoLC solution): %.3f%%\n', 100*stats.FractionInvalid);
    fprintf('  detections      : %d (%.4f%%)\n', ...
        stats.NumDetections, 100*stats.NumDetections/numel(I));
    fprintf('  log vs amplitude domain mismatch: %d px (%.6f%%)\n', ...
        logdomain.NumMismatch, logdomain.PercentMismatch);
end
end


function prm = apply_rmax(prm, rMax)
%APPLY_RMAX  Optional extra upper bound on r, for reproducing the legacy
%   `3*c3^2 <= 8*c2^3` (r <= 8/3) guard exactly.
if isfinite(rMax)
    prm.valid = prm.valid & (prm.r <= rMax);
    prm.FractionInvalid = 1 - mean(prm.valid(:));
    prm.RMax = rMax;
end
end
