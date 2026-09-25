function [detection_map, threshold_map, prm, stats, logdomain] = ...
    KCFAR_Floating(I, sli, guard, Pfa, varargin)
%KCFAR_FLOATING  Single-pass windowed K-distribution CFAR detector.
%
%   Floating-point reference model for the streaming FPGA implementation,
%   same contract as WeibullCFAR_Floating.m / BurrCFAR_Floating.m. Composed
%   from
%       cfar_front_end  -- shared with all six detectors
%       KCFAR_Params    -- c2 -> shape `a` (single-look, L=1) via the
%                          inverse-trigamma MoLC relation
%       KCFAR_TLog      -- delta = [quantile of the K-amplitude survival
%                          mixture integral] - kappa1(a)
%
%   NOT the same pipeline as CFAR K/legacy/Main_Estimation_and_Detection_K.m
%   -----------------------------------------------------------------------
%   That legacy script runs Gamma-MAP texture filtering FIRST, then a
%   separate Gamma-texture CFAR stage (nkgmolc.m) on the filtered image --
%   a two-stage pipeline with no equivalent in this project's shared
%   single-pass front end, and not a fair comparison against the other five
%   detectors (which all see raw single-look log-amplitude data with no
%   despeckling stage). This file instead applies the true K-distribution
%   MoLC/quantile relations directly to the same cfar_front_end moments
%   every other detector uses, so the six-way comparison stays apples-to-
%   apples.
%
%   NAME-VALUE OPTIONS
%     'AMin' (0.05), 'AMax' (100) : shape-parameter clamps, see KCFAR_Params
%     'FrontEnd', 'FrontEndOpts', 'Verbose'
%
%   See also: cfar_front_end, cfar_compose, KCFAR_Params, KCFAR_TLog

p = inputParser;
addParameter(p, 'AMin', 0.05);
addParameter(p, 'AMax', 100);
addParameter(p, 'FrontEnd',     [],  @(v) isempty(v) || isstruct(v));
addParameter(p, 'FrontEndOpts', {},  @iscell);
addParameter(p, 'Verbose',      false, @(v) islogical(v) && isscalar(v));
parse(p, varargin{:});

if isempty(p.Results.FrontEnd)
    feArg = p.Results.FrontEndOpts;
else
    feArg = p.Results.FrontEnd;
end

paramsFcn = @(c2, c3, skew) KCFAR_Params(c2, 'AMin', p.Results.AMin, 'AMax', p.Results.AMax);
tlogFcn   = @(prm, pf) KCFAR_TLog(prm, pf);

[detection_map, threshold_map, prm, stats, logdomain] = ...
    cfar_compose(I, sli, guard, Pfa, paramsFcn, tlogFcn, 'K', feArg);

if p.Results.Verbose
    v = prm.valid;
    fprintf('KCFAR_Floating: %dx%d, sli=%d guard=%d N=%d Pfa=%g\n', ...
        size(I,1), size(I,2), sli, guard, stats.NumReferenceCells, Pfa);
    if any(v(:))
        fprintf('  a (shape) : min=%.4f median=%.4f max=%.4f\n', ...
            min(prm.a(v)), median(prm.a(v)), max(prm.a(v)));
    end
    fprintf('  INVALID windows (no K MoLC solution, c2 <= pi^2/24): %.3f%%\n', ...
        100*stats.FractionInvalid);
    fprintf('  detections : %d (%.4f%%)\n', ...
        stats.NumDetections, 100*stats.NumDetections/numel(I));
    fprintf('  log vs amplitude domain mismatch: %d px (%.6f%%)\n', ...
        logdomain.NumMismatch, logdomain.PercentMismatch);
end
end
