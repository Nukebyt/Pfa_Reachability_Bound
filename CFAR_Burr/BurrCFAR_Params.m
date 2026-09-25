function prm = BurrCFAR_Params(c2, c3, varargin)
%BURRCFAR_PARAMS  Method-of-log-cumulants estimator for Burr-XII clutter.
%
%   prm = BURRCFAR_PARAMS(c2, c3, ...)
%
%   ---------------------------------------------------------------------
%   THE MoLC SYSTEM
%   ---------------------------------------------------------------------
%   From CFAR_Burr/legacy/burrfun.m and burrmolc.m, with rho the first shape
%   parameter, kappa the second and eta the scale:
%       c2 = ( psi(1,kappa) + psi(1,1) ) / rho^2
%       c3 = ( psi(2,1) - psi(2,kappa) ) / rho^3
%       eta = exp( c1 - (psi(1) - psi(kappa)) / rho )
%
%   ---------------------------------------------------------------------
%   THE 2x2 SYSTEM IS EXACTLY 1-D -- eliminate rho analytically
%   ---------------------------------------------------------------------
%   From the first equation  rho^2 = (psi(1,kappa)+psi(1,1))/c2, and from
%   the second  rho^3 = (psi(2,1)-psi(2,kappa))/c3. Dividing the second by
%   the first and equating to rho = sqrt(...) removes rho entirely:
%
%       psi(2,1) - psi(2,kappa)
%       -----------------------------------  =  c3 / c2^(3/2)  =:  s
%       ( psi(1,kappa) + psi(1,1) )^(3/2)
%
%   The right-hand side is the DIMENSIONLESS log-skewness -- scale-free,
%   exactly as for the generalized gamma's r = c3^2/c2^3. So kappa depends on
%   ONE scalar, and the ostensibly 2-D (c2,c3) shape table collapses to a 1-D
%   ROM in Phase 3. Given kappa, rho follows from c2 with one square root.
%
%   Call the left-hand side B(kappa). Measured over kappa in [1e-4, 1e4],
%   B is strictly decreasing with range
%
%       B(kappa) -> 2         as kappa -> 0+
%       B(1)      = 0                             (psi(2,1) cancels)
%       B(kappa) -> -1.139443 as kappa -> inf
%
%   so a Burr-XII MoLC solution exists if and only if
%
%           -1.139443 < s < 2                      <-- the support condition
%
%   Windows outside that band have no solution and are counted in
%   FractionInvalid rather than being handed to a solver. (The legacy code
%   called fsolve from a fixed t0 = [0.8 0.8] with no support test at all,
%   so an out-of-support window simply produced whatever fsolve stopped at,
%   flagged only by an exit code the main script recorded but never acted on.)
%
%   Because B is monotone, the inversion is done with mono_table_invert --
%   a dense monotone table plus one secant polish. That is both a replacement
%   for the unavailable fsolve and the floating-point model of the Phase 3
%   ROM; they differ only in table density and output word length.
%
%   ---------------------------------------------------------------------
%   RHO, AND WHY ITS SIGN IS NOT AMBIGUOUS
%   ---------------------------------------------------------------------
%       rho = sqrt( (psi(1,kappa) + psi(1,1)) / c2 )
%   psi(1,.) > 0 everywhere on (0,inf) and c2 > 0, so the argument is always
%   positive and rho is real. The positive root is the correct one: rho is a
%   Burr-XII shape parameter and must be > 0, and the sign information in c3
%   has already been consumed by the sign of s when solving for kappa.
%
%   NAME-VALUE OPTIONS
%     'KappaMin' (1e-3), 'KappaMax' (1e3) : clamps on kappa
%     'RhoMin'   (0.05), 'RhoMax'  (50)   : clamps on rho
%
%   OUTPUT (struct prm)
%     kappa, rho        : HxW clamped parameter maps
%     kappa_raw, rho_raw : unclamped (diagnostic / Phase 3 LUT range sizing)
%     s                 : HxW c3/c2^1.5, the 1-D LUT address variable
%     valid             : HxW logical, from the support condition above
%     FractionInvalid, FractionAtMin, FractionAtMax
%
%   See also: mono_table_invert, BurrCFAR_TLog, BurrCFAR_Floating, cfar_front_end

p = inputParser;
addParameter(p, 'KappaMin', 1e-3, @(v) isnumeric(v) && isscalar(v) && v > 0);
addParameter(p, 'KappaMax', 1e3,  @(v) isnumeric(v) && isscalar(v) && v > 0);
addParameter(p, 'RhoMin',   0.05, @(v) isnumeric(v) && isscalar(v) && v > 0);
addParameter(p, 'RhoMax',   50,   @(v) isnumeric(v) && isscalar(v) && v > 0);
parse(p, varargin{:});
kMin = p.Results.KappaMin;
kMax = p.Results.KappaMax;
rMin = p.Results.RhoMin;
rMax = p.Results.RhoMax;

%% ---- The scale-free address variable -----------------------------------
s = c3 ./ c2.^1.5;

%% ---- Invert B(kappa) = s -----------------------------------------------
persistent kgrid
if isempty(kgrid)
    kgrid = logspace(-4, 4, 40001)';
end

PSI11 = psi(1,1);
PSI21 = psi(2,1);
Bfun  = @(kap) (PSI21 - psi(2, kap)) ./ (psi(1, kap) + PSI11).^1.5;

[kappa_raw, valid] = mono_table_invert(Bfun, kgrid, s);

kappa = min(max(kappa_raw, kMin), kMax);

%% ---- rho from kappa and c2 ---------------------------------------------
rho_raw = sqrt( (psi(1, kappa) + PSI11) ./ c2 );
rho     = min(max(rho_raw, rMin), rMax);

valid = valid & isfinite(kappa_raw) & isfinite(rho_raw) & rho_raw > 0;

kappa(~valid) = NaN;
rho(~valid)   = NaN;

%% ---- Output ------------------------------------------------------------
prm = struct();
prm.Name            = 'BurrXII';
prm.kappa           = kappa;
prm.kappa_raw       = kappa_raw;
prm.rho             = rho;
prm.rho_raw         = rho_raw;
prm.s               = s;
prm.valid           = valid;
prm.FractionInvalid = 1 - mean(valid(:));
if any(valid(:))
    prm.FractionAtMin = mean(kappa_raw(valid) <= kMin + 1e-12);
    prm.FractionAtMax = mean(kappa_raw(valid) >= kMax - 1e-12);
else
    prm.FractionAtMin = NaN;
    prm.FractionAtMax = NaN;
end
prm.KappaMin = kMin;  prm.KappaMax = kMax;
prm.RhoMin   = rMin;  prm.RhoMax   = rMax;
end
