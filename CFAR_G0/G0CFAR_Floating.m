function [detection_map, threshold_map, prm, stats, logdomain] = ...
    G0CFAR_Floating(I, sli, guard, Pfa, varargin)
%G0CFAR_FLOATING  Single-pass windowed G0 CFAR detector.
%
%   Floating-point reference model for the streaming FPGA implementation,
%   with the same contract as WeibullCFAR_Floating.m. Composed from
%       cfar_front_end   -- shared with all five detectors
%       G0CFAR_Params    -- (c2,c3) -> L, alpha  by bounded 1-D bisection
%       G0CFAR_TLog      -- delta = 0.5*[psi(u) - psi(L) + log(x/(1-x))]
%
%   This is one of the two detectors slated for FPGA implementation first.
%
%   DIFFERENCES FROM CFAR_G0/legacy/Main_CFAR_G0.m
%   ----------------------------------------------
%   1. Vectorised single pass instead of the per-pixel double loop.
%   2. fsolve (no Optimization Toolbox installed) is replaced by an exact
%      reduction of the 2x2 MoLC system to one bounded, monotone, scalar
%      residual plus a closed-form support test -- see G0CFAR_Params.m. This
%      removes the dependence on fsolve's fixed starting guess and makes
%      "no solution exists" a detected condition rather than a silent
%      convergence failure.
%   3. finv (also Statistics Toolbox) is replaced by the exact incomplete-beta
%      identity -- see G0CFAR_TLog.m.
%   4. The legacy CA-CFAR normalisation (`condi = mean + std*thres`) is
%      dropped, uniformly across all five detectors, so Phase 2 compares
%      like with like. See LognormalCFAR_Floating.m for the reasoning.
%   5. The legacy script pads with the constant 10 (`padarray(I,...,10,...)`)
%      where the other legacy scripts pad with max(I(:)); cfar_front_end uses
%      symmetric padding for all five, which is both convention-independent
%      and closest to what a line-buffer boundary policy can do.
%   6. legacy/thres_G0.m + thres_G0_fun.m computed the threshold by calling
%      fsolve on the G0 CDF written with a generalized hypergeometric
%      function (hypergeomq.m, a Symbolic-Math-dependent helper). That path
%      is not used here: the incomplete-beta form is exact, closed form, and
%      needs no symbolic toolbox. hypergeomq.m is retained in legacy/ only
%      as a reference for that alternative derivation.
%
%   NAME-VALUE OPTIONS
%     'Mode'  : 'LA' (default, estimates both L and alpha, matching
%               legacy/g0molc.m) | 'L1' (single-look, L fixed -- the cheap,
%               1-D-ROM hardware variant; expect a high invalid fraction,
%               see G0CFAR_Params.m)
%     'L'     : fixed look count for 'L1' (default 1)
%     'Iters' : bisection iterations for 'LA' (default 30)
%     'UMin' (0.05), 'UMax' (200), 'LMax' (100)
%     'FrontEnd', 'FrontEndOpts', 'Verbose'
%
%   See also: cfar_front_end, cfar_compose, G0CFAR_Params, G0CFAR_TLog

p = inputParser;
addParameter(p, 'Mode',         'LA');
addParameter(p, 'L',            1);
addParameter(p, 'Iters',        30);
addParameter(p, 'UMin',         0.05);
addParameter(p, 'UMax',         200);
addParameter(p, 'LMax',         100);
addParameter(p, 'FrontEnd',     [],  @(v) isempty(v) || isstruct(v));
addParameter(p, 'FrontEndOpts', {},  @iscell);
addParameter(p, 'Verbose',      false, @(v) islogical(v) && isscalar(v));
parse(p, varargin{:});

if isempty(p.Results.FrontEnd)
    feArg = p.Results.FrontEndOpts;
else
    feArg = p.Results.FrontEnd;
end

paramsFcn = @(c2, c3, skew) G0CFAR_Params(c2, c3, ...
    'Mode',  p.Results.Mode,  'L',    p.Results.L, ...
    'Iters', p.Results.Iters, 'UMin', p.Results.UMin, ...
    'UMax',  p.Results.UMax,  'LMax', p.Results.LMax);
tlogFcn = @(prm, pf) G0CFAR_TLog(prm, pf);

[detection_map, threshold_map, prm, stats, logdomain] = ...
    cfar_compose(I, sli, guard, Pfa, paramsFcn, tlogFcn, 'G0', feArg);

if p.Results.Verbose
    v = prm.valid;
    fprintf('G0CFAR_Floating: %dx%d, sli=%d guard=%d N=%d Pfa=%g mode=%s\n', ...
        size(I,1), size(I,2), sli, guard, stats.NumReferenceCells, Pfa, prm.Mode);
    if any(v(:))
        fprintf('  u = -alpha : min=%.3f median=%.3f max=%.3f\n', ...
            min(prm.u(v)), median(prm.u(v)), max(prm.u(v)));
        fprintf('  L          : min=%.3f median=%.3f max=%.3f\n', ...
            min(prm.L(v)), median(prm.L(v)), max(prm.L(v)));
    end
    fprintf('  INVALID windows (no G0 MoLC solution): %.3f%%\n', 100*stats.FractionInvalid);
    fprintf('  detections : %d (%.4f%%)\n', ...
        stats.NumDetections, 100*stats.NumDetections/numel(I));
    fprintf('  log vs amplitude domain mismatch: %d px (%.6f%%)\n', ...
        logdomain.NumMismatch, logdomain.PercentMismatch);
end
end
