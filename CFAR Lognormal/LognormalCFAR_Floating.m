function [detection_map, threshold_map, prm, stats, logdomain] = ...
    LognormalCFAR_Floating(I, sli, guard, Pfa, varargin)
%LOGNORMALCFAR_FLOATING  Single-pass windowed Lognormal CFAR detector.
%
%   [detection_map, threshold_map, prm, stats, logdomain] = ...
%       LognormalCFAR_Floating(I, sli, guard, Pfa, ...)
%
%   Floating-point reference model for the streaming FPGA implementation,
%   with the same contract as WeibullCFAR_Floating.m. Built from three
%   pieces:
%       cfar_front_end          -- shared with all five detectors
%       LognormalCFAR_Params    -- mu = c1, sigma = sqrt(c2)
%       LognormalCFAR_TLog      -- delta = z(Pfa) * sigma
%   composed by cfar_compose, which also runs the amplitude-vs-log-domain
%   equivalence check on every call.
%
%   DIFFERENCES FROM CFAR Lognormal/legacy/Main_CFAR_Lognormal.m
%   -----------------------------------------------------------
%   1. Single pass, fully vectorised (integral-image window sums) instead of
%      a per-pixel double loop with nonzeros() -- the original is O(h*w*sli^2)
%      and takes minutes per image at the window sizes Phase 2 sweeps.
%   2. The original applies a CA-CFAR normalisation on top of the lognormal
%      threshold:
%            condi = mean(clutter) + std(clutter)*thres;   I > condi
%      That double-counts the clutter level -- `thres` is already an absolute
%      amplitude threshold derived from the fitted distribution, so
%      multiplying it by the clutter std and adding the clutter mean destroys
%      the calibration that makes the nominal Pfa mean anything. Every model
%      in this project uses the direct comparison instead, uniformly, so the
%      Phase 2 cross-detector numbers are comparable. The original form is
%      preserved verbatim in legacy/ if a side-by-side is ever wanted.
%   3. The original's guard-hole indexing (`temp(cen-t:cen+t, ...)` with
%      t = tsli-tguard) makes the excluded region depend on the DIFFERENCE of
%      the half-widths rather than on the guard size itself, so the actual
%      hole is (sli-guard+1) wide, not guard wide. cfar_front_end uses the
%      conventional definition (hole = guard x guard, N = sli^2 - guard^2)
%      which is also what the Weibull RTL implements.
%
%   INPUTS
%     I     : HxW SAR intensity image, 8-bit-valued
%     sli   : sliding window size (odd, >= 3)
%     guard : guard region size (odd, >= 1, < sli)
%     Pfa   : desired probability of false alarm, 0 < Pfa < 0.5
%
%   NAME-VALUE OPTIONS
%     'SigmaMin' (0.05), 'SigmaMax' (2.0)  -- clamps, see LognormalCFAR_Params
%     'FrontEnd' -- a precomputed cfar_front_end struct to reuse
%     'FrontEndOpts' -- cell array forwarded to cfar_front_end
%     'Verbose'  -- print a summary (default false)
%
%   See also: cfar_front_end, cfar_compose, LognormalCFAR_Params,
%             LognormalCFAR_TLog

p = inputParser;
addParameter(p, 'SigmaMin',     0.05, @(v) isnumeric(v) && isscalar(v) && v > 0);
addParameter(p, 'SigmaMax',     2.0,  @(v) isnumeric(v) && isscalar(v) && v > 0);
addParameter(p, 'FrontEnd',     [],   @(v) isempty(v) || isstruct(v));
addParameter(p, 'FrontEndOpts', {},   @iscell);
addParameter(p, 'Verbose',      false, @(v) islogical(v) && isscalar(v));
parse(p, varargin{:});

if isempty(p.Results.FrontEnd)
    feArg = p.Results.FrontEndOpts;
else
    feArg = p.Results.FrontEnd;
end

paramsFcn = @(c2, c3, skew) LognormalCFAR_Params(c2, ...
    'SigmaMin', p.Results.SigmaMin, 'SigmaMax', p.Results.SigmaMax);
tlogFcn   = @(prm, pf) LognormalCFAR_TLog(prm, pf);

[detection_map, threshold_map, prm, stats, logdomain] = ...
    cfar_compose(I, sli, guard, Pfa, paramsFcn, tlogFcn, 'Lognormal', feArg);

if p.Results.Verbose
    fprintf('LognormalCFAR_Floating: %dx%d, sli=%d guard=%d N=%d Pfa=%g\n', ...
        size(I,1), size(I,2), sli, guard, stats.NumReferenceCells, Pfa);
    fprintf('  sigma (clamped) : min=%.4f mean=%.4f max=%.4f  [%.2f, %.2f]\n', ...
        min(prm.sigma(:)), mean(prm.sigma(:)), max(prm.sigma(:)), ...
        prm.SigmaMin, prm.SigmaMax);
    fprintf('  clamp saturation: %.3f%% low, %.3f%% high | invalid windows: %.3f%%\n', ...
        100*prm.FractionAtMin, 100*prm.FractionAtMax, 100*stats.FractionInvalid);
    fprintf('  detections      : %d (%.4f%%)\n', ...
        stats.NumDetections, 100*stats.NumDetections/numel(I));
    fprintf('  log vs amplitude domain mismatch: %d px (%.6f%%)\n', ...
        logdomain.NumMismatch, logdomain.PercentMismatch);
end
end
