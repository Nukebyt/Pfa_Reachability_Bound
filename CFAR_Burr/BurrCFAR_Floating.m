function [detection_map, threshold_map, prm, stats, logdomain] = ...
    BurrCFAR_Floating(I, sli, guard, Pfa, varargin)
%BURRCFAR_FLOATING  Single-pass windowed Burr-XII CFAR detector.
%
%   Floating-point reference model for the streaming FPGA implementation,
%   with the same contract as WeibullCFAR_Floating.m. Composed from
%       cfar_front_end    -- shared with all five detectors
%       BurrCFAR_Params   -- (c2,c3) -> kappa, rho via the scale-free
%                            skewness s = c3/c2^1.5
%       BurrCFAR_TLog     -- delta = [psi(kappa) - psi(1) + log(Pfa^-1/kappa - 1)]/rho
%
%   DIFFERENCES FROM CFAR_Burr/legacy/Main_CFAR_BurrXII.m
%   -----------------------------------------------------
%   1. Vectorised single pass instead of the per-pixel double loop.
%   2. fsolve on the 2x2 MoLC system (no Optimization Toolbox installed) is
%      replaced by an EXACT analytic reduction to one monotone scalar
%      equation in kappa plus a closed-form support condition
%      (-1.139443 < s < 2) -- see BurrCFAR_Params.m. The legacy code solved
%      from a fixed t0 = [0.8 0.8] with no support test; its exit flag was
%      stored in flag_check but never acted on, so out-of-support windows
%      contributed whatever fsolve happened to stop at.
%   3. The legacy script's threshold comparison is `I_pad(i,j) > thres`, i.e.
%      it computes `condi = meu + sigma*thres` on the line above and then
%      never uses it -- the one legacy script of the four that does NOT
%      apply the CA-CFAR normalisation. All five models here use the direct
%      comparison, so Burr is the detector whose legacy decision rule this
%      model is closest to.
%   4. The legacy guard-hole indexing is `temp(cen-tguard:cen+tguard,...)`,
%      which excludes a guard x guard hole -- matching cfar_front_end's
%      convention (unlike the Lognormal and G0 legacy scripts, which use
%      cen-t:cen+t with t = tsli-tguard and so exclude a hole of a different
%      size than `guard`).
%   5. Pfa^(-1/kappa) is never formed directly; see the overflow discussion
%      in BurrCFAR_TLog.m. At Pfa = 1e-6 the legacy expression overflows to
%      Inf for kappa below ~0.02 and silently suppresses all detection there.
%
%   NAME-VALUE OPTIONS
%     'KappaMin' (1e-3), 'KappaMax' (1e3), 'RhoMin' (0.05), 'RhoMax' (50)
%     'FrontEnd', 'FrontEndOpts', 'Verbose'
%
%   See also: cfar_front_end, cfar_compose, BurrCFAR_Params, BurrCFAR_TLog

p = inputParser;
addParameter(p, 'KappaMin',     1e-3);
addParameter(p, 'KappaMax',     1e3);
addParameter(p, 'RhoMin',       0.05);
addParameter(p, 'RhoMax',       50);
addParameter(p, 'FrontEnd',     [],  @(v) isempty(v) || isstruct(v));
addParameter(p, 'FrontEndOpts', {},  @iscell);
addParameter(p, 'Verbose',      false, @(v) islogical(v) && isscalar(v));
parse(p, varargin{:});

if isempty(p.Results.FrontEnd)
    feArg = p.Results.FrontEndOpts;
else
    feArg = p.Results.FrontEnd;
end

paramsFcn = @(c2, c3, skew) BurrCFAR_Params(c2, c3, ...
    'KappaMin', p.Results.KappaMin, 'KappaMax', p.Results.KappaMax, ...
    'RhoMin',   p.Results.RhoMin,   'RhoMax',   p.Results.RhoMax);
tlogFcn = @(prm, pf) BurrCFAR_TLog(prm, pf);

[detection_map, threshold_map, prm, stats, logdomain] = ...
    cfar_compose(I, sli, guard, Pfa, paramsFcn, tlogFcn, 'BurrXII', feArg);

if p.Results.Verbose
    v = prm.valid;
    fprintf('BurrCFAR_Floating: %dx%d, sli=%d guard=%d N=%d Pfa=%g\n', ...
        size(I,1), size(I,2), sli, guard, stats.NumReferenceCells, Pfa);
    if any(v(:))
        fprintf('  kappa : min=%.4f median=%.4f max=%.4f\n', ...
            min(prm.kappa(v)), median(prm.kappa(v)), max(prm.kappa(v)));
        fprintf('  rho   : min=%.4f median=%.4f max=%.4f\n', ...
            min(prm.rho(v)), median(prm.rho(v)), max(prm.rho(v)));
        fprintf('  s = c3/c2^1.5 : min=%.4f median=%.4f max=%.4f  [support: -1.1394 .. 2]\n', ...
            min(prm.s(v)), median(prm.s(v)), max(prm.s(v)));
    end
    fprintf('  INVALID windows (no Burr MoLC solution): %.3f%%\n', 100*stats.FractionInvalid);
    fprintf('  detections : %d (%.4f%%)\n', ...
        stats.NumDetections, 100*stats.NumDetections/numel(I));
    fprintf('  log vs amplitude domain mismatch: %d px (%.6f%%)\n', ...
        logdomain.NumMismatch, logdomain.PercentMismatch);
end
end
