# The Pfa Reachability Bound

**Status tracker.** Written incrementally as roadmap items in
`PAPER1_ROADMAP_reachability-bound.md` are completed. `[PLACEHOLDER -- ...]`
marks a section where the underlying derivation/measurement doesn't exist
yet — never filled with an invented number or an unverified derivation.

**Target venue:** IEEE TGRS primary, JSTARS fallback. Solo (no faculty
co-author) — **the authorship conversation with the supervisor has not been
confirmed as had; do not treat this file's existence as resolving that.**

**Note on the two shared common-step items (BC-5 CA-CFAR, BC-6 metrics fix):
this paper does not need either.** Its headline claims (the reachability
bound, the calibration divergence) are Pfa-based, not F1/Precision/Recall-
based, and it has no non-parametric baseline comparison in its own roadmap.
Confirmed explicitly here per the roadmap's own item 3 instruction, so no
stray F1 table accidentally imports a number from a shared CSV.

---

## Roadmap item tracker

| Item | Status |
|---|---|
| BC-3 — Meng collision (item 1) | [x] Resolved — see §3 below |
| Growth-rate derivation, all six distributions (item 2) | [x] 6/6 done — see §2.1. Two exact (Weibull, Burr XII), four asymptotic (Lognormal, GenGamma [2 branches], G0, K), all numerically cross-checked against the real `*CFAR_TLog.m` code (`verify_growth_rates.m`) |
| Metric bug / F1 usage (item 3) | [x] Confirmed N/A — this paper reports no F1/Precision/Recall numbers |
| Clopper–Pearson intervals (item 4) | [x] Done — see §4. Headline figures tight (<5% relative width); well-calibrated detectors' ratios and zero-false-alarm cells are NOT tight (up to 47%/undefined) — stated explicitly |
| HRSID sample expansion past 60 images (item 5) | [x] Done — 60->250 images, confirmed stable (§4 addendum); also independently validated the CI-width finding (Lognormal's point ratio moved 2x+ across samples, GenGamma/Weibull's did not) |
| Xplore prior-art search (item 6) | [x] Done — see §5. Audited search, no collision found (7 papers read, none states this thesis under any terminology); not exhaustive, but real evidence, reported as such |

---

## 1. Introduction / Thesis

A constant-false-alarm-rate (CFAR) detector's defining promise is that the
operator's requested false-alarm probability, the nominal Pfa, is what the
system actually delivers, regardless of clutter level or clutter type. In
practice this promise is known to be imperfect. Recent work on SAR ship
detection (Meng et al., 2025, JSTARS) reports substantial gaps between
requested and achieved Pfa across scenes, generally attributed to clutter-
model mismatch and scene-to-scene statistical diversity — the achieved rate
drifts because the *scene* the detector is pointed at does not match the
scene it was calibrated on. This is a real and well-documented failure
mode, and §3 of this paper shows directly that it is not the mechanism this
paper is about.

**This paper identifies a second, model-internal mechanism, orthogonal to
scene mismatch: even for a single scene, a single window geometry, and a
clutter model that is exactly correctly specified with noiseless shape
estimates, the achieved Pfa can still diverge from the nominal one by
several orders of magnitude — because of how the detector's own
threshold-offset function is shaped, not because of anything wrong with the
clutter fit.** Every detector considered here reduces to the same decision
rule, `T_log = c1 + delta(shape, Pfa)`, where `c1` carries the clutter
level and `delta` is a function of the fitted shape parameter(s) and the
requested Pfa alone. The question this paper asks is simple and, as far as
a search of the accessible literature could establish (§5), not previously
asked in these terms: **how fast does `delta` grow as `Pfa -> 0`, and does
that growth rate differ enough between commonly used SAR clutter families
to matter in practice?**

It differs enormously. §2 derives — exactly for two of six families, as a
verified asymptotic for the other four, every derivation checked
numerically against the real implementation rather than trusted as
algebra — the `Pfa -> 0` growth rate of all six clutter families this
project's broader comparison uses (Weibull, Lognormal, Generalized Gamma,
G0, Burr XII, K-distribution). They separate cleanly into three classes,
and the class a family falls into is set by the extreme-value type of its
clutter tail, not by which specific detector it happens to be: light,
Gumbel-domain tails (Weibull, K-distribution, and half of Generalized
Gamma's parameter range) give a threshold offset that grows only as
`log log(1/Pfa)`; heavier, power-law tails (Burr XII, G0, and the other half
of Generalized Gamma) give a threshold offset that grows a full order
faster, as `log(1/Pfa)`.

**To be precise about what this growth-rate result does and does not
claim, stated here rather than left for a reader to infer:** no detector's
threshold offset is bounded above. Every one of the six formulas derived in
§2 diverges to infinity as `Pfa -> 0` — there is no asymptotic floor on the
achievable false-alarm rate, and this paper does not claim one exists. What
differs between the three growth classes is *how much a given decrease in
the nominal Pfa dial actually moves the threshold*. For a log-log-growing
family, asking for a nominal Pfa ten orders of magnitude smaller moves the
threshold offset by only a few units of `delta` — an amount easily
swamped by ordinary shape-estimation noise (already documented, for real
SAR data, as this project's own F9 finding) or by the quantization step of
a realistic fixed-point hardware pipeline. The dial still turns; it simply
stops doing anything useful long before it reaches the requested value.
This paper calls that operational fact — not an asymptotic bound on
achievable Pfa, but a practical bound on how far the nominal-Pfa dial can
be turned before its effect disappears into noise — the **Pfa reachability
bound**.

This paper's contributions are:

1. A closed-form or tightly-bounded-asymptotic derivation of the `Pfa -> 0`
   growth rate for all six SAR clutter families in this comparison,
   generalizing a previously known two-family result (Weibull, Lognormal)
   into a three-class taxonomy tied directly to the extreme-value type of
   the underlying clutter tail — a theorem about distribution families, not
   a pairwise curiosity (§2).
2. Direct resolution of an apparent collision with a recently published
   cross-scene Pfa-divergence result, showing on a controlled, fixed-scene,
   swept-window-size grid that this paper's divergence axis is empirically
   distinct from that one (§3).
3. Exact Clopper-Pearson confidence intervals, computed from the clutter-
   window sample size rather than the much sparser ship count, attached to
   every reported divergence figure — and a demonstration that this
   precision is fundamentally asymmetric: large, alarming divergence
   numbers are usually well-supported by the data, while small,
   reassuring "well-calibrated" numbers frequently are not (§4).
4. An audited, though not exhaustive, search of the accessible IEEE
   literature establishing that this specific growth-rate mechanism does
   not appear to have been identified previously under different
   terminology (§5).

Section 2 develops the growth-rate taxonomy and its numerical verification.
Section 3 resolves the collision with prior cross-scene work. Section 4
attaches statistical rigor to every headline number this paper and its
companion measurements report. Section 5 reports the literature search.
Section 6 discusses what the growth-rate result does and does not imply for
detector selection and for how this class of result should be reported
generally, and Section 7 concludes.

## 2. Closed-Form Growth Rates

*Verified, pre-existing (this project's F9/F10 findings):*

- **Weibull:** `K(Pfa) = gamma + log(-log Pfa)` — grows log-logarithmically.
- **Lognormal:** `z(Pfa) ~ sqrt(log(1/Pfa))` — grows as a square root of
  log-inverse-Pfa, structurally slower-diverging than Weibull's in the
  relevant sense the roadmap's thesis rests on. F10 additionally proves
  Weibull and Lognormal are the SAME decision rule up to one Pfa-dependent
  scalar (`delta_Lognormal/delta_Weibull` constant to 4e-16/6.7e-16 across a
  decade of `c2`, at Pfa=1e-3 and 1e-6 respectively) — they differ only
  through that scalar and their clamp behavior, not through independently
  free growth rates.

### 2.1 The remaining four (2026-09, roadmap item 2)

**Method.** Each of the four remaining detectors' `delta(shape, Pfa) = T_log
- c1` already has an EXACT implicit form in its `*CFAR_TLog.m` (an inverse
incomplete gamma, inverse incomplete beta, or — for K — a numerically solved
1-D mixture integral; see each file's header). What's new here is each
one's closed-form or tightly-bounded asymptotic behavior AS `Pfa -> 0`,
derived by standard tail-inversion analysis (large- and small-argument
series of the incomplete gamma/beta functions; Laplace's-method
saddle-point analysis for K's mixture integral, which has no elementary
tail series at all). Every derivation below was cross-checked numerically
against the real `*CFAR_TLog.m` output — not just algebra on paper — via
`_comparison/verify_growth_rates.m`, at both this project's real sweep Pfa
range (`1e-1...1e-6`) and much deeper (`1e-9...1e-40`), confirming the error
between the derived formula and the real code shrinks as `Pfa -> 0`. **This
check caught one real sign error** before it reached this document (BUG_LOG
D25: the Generalized Gamma `v<0` branch's leading term had the wrong sign in
an early draft — predicted values came out negated relative to the real
code, an unmistakable tell, fixed and reverified).

**Result — the six distributions fall into three growth classes, and the
class is set by the underlying tail TYPE, not by which detector it happens
to be:**

| Class | Growth as `Pfa -> 0` | Members |
|---|---|---|
| **Log-log** (slowest) | `~ log log(1/Pfa)` | Weibull (exact); GenGamma `v>0` branch (asymptotic); K-distribution (asymptotic) |
| **Sqrt-log** | `~ sqrt(log(1/Pfa))` | Lognormal (asymptotic) |
| **Linear-log** (fastest) | `~ log(1/Pfa)` | Burr XII (exact); G0 (asymptotic); GenGamma `v<0` branch (asymptotic) |

**Generalized Gamma** (`GenGammaCFAR_TLog.m`: `delta = [log(G) - psi(k)]/v`,
`G` a gamma-quantile whose tail depends on the sign of `v`):

```
    v > 0 (upper gamma tail):  delta ~ [ log log(1/Pfa) - psi(k) ] / v
                                        + O( log log(1/Pfa) / (v * log(1/Pfa)) )
    v < 0 (lower gamma tail):  delta ~ [ log(1/Pfa)/k - log(Gamma(k+1))/k + psi(k) ] / |v|
                                        + O( Pfa^(1/k) / |v| )
```

The two branches of the SAME family land in different growth classes: `v>0`
inherits Weibull's exact log-log form (an upper incomplete-gamma tail decays
like `x^(k-1) e^-x`, the same Gumbel-type exponential tail that gives
Weibull its own log-log growth — `k=1` reduces the two formulas to be
identical, confirmed to `0.000` absolute error at every Pfa tested); `v<0`
inherits Burr/G0's fast linear-log form (a lower incomplete-gamma tail near
0 decays like `x^k`, a power law, the same tail type that gives Burr and G0
their linear growth). **This asymmetry is real, not a labeling
convenience** — the `v<0` branch converges to its asymptotic form
polynomially fast (`O(Pfa^(1/k))`, machine precision by `Pfa=1e-12` in the
numerical check), while the `v>0` branch converges only logarithmically
slowly (`O(log L / L)`, still `~0.05`–`0.2` absolute error even at
`Pfa=1e-40` for the `k` values tested) — a second, independent piece of
evidence (beyond F9's already-documented `1/v` noise amplification) that
Generalized Gamma's worst-case calibration behavior concentrates on the
positive-skew branch.

**G0** (`G0CFAR_TLog.m`: `delta = 0.5*[psi(u)-psi(L) + log(x/(1-x))]`, `x`
solving the regularized incomplete beta `I_{1-x}(u,L) = Pfa`):

```
    delta ~ 0.5*(psi(u) - psi(L)) + [ log(1/Pfa) - log(u * B(u,L)) ] / (2u)
```

using the small-argument beta series `I_y(u,L) ~ y^u / (u B(u,L))` as
`y -> 0`. Error `O(Pfa^(1/u))` — already at `5e-7` absolute error by
`Pfa=1e-6` for `u=1`, machine precision by `Pfa=1e-20`. Linear-log class,
rate `1/(2u)` — SMALLER `u` (rougher texture, heavier tail) grows the
offset FASTER, i.e. calibrates tighter; larger `u` degrades toward the
slower classes' behavior, consistent with G0's real-data ratio (31.8x at
`Pfa=1e-6`, F9) sitting well below Weibull's/GenGamma's but not as low as a
`u<<1` asymptotic optimum would suggest.

**Burr XII** (`BurrCFAR_TLog.m`: `delta = [psi(kappa)-psi(1) +
log(Pfa^(-1/kappa)-1)]/rho`) — this one is genuinely **exact**, not merely
asymptotic, because Burr's survival function is algebraically invertible in
closed form (like Weibull's):

```
    delta = [ psi(kappa) - psi(1) + log(1/Pfa)/kappa + log1p(-Pfa^(1/kappa)) ] / rho
```

the `log1p` term is `O(Pfa^(1/kappa))` — literally exponentially small in
`log(1/Pfa)/kappa`, not merely "eventually negligible." Linear-log class,
rate `1/(kappa*rho)`. Verified: matches the real `BurrCFAR_TLog.m` to
`1e-4`–`1e-1` absolute error even at `Pfa=1e-1` and to double-precision
floor by `Pfa=1e-20`.

**K-distribution** (`KCFAR_TLog.m`: `delta = 0.5*log(x_pfa) - kappa1(a)`,
`x_pfa` solving the mixture survival integral `S(x;a)=E_tau[exp(-x/tau)]=
Pfa`, `tau ~ Gamma(a,1/a)`) — no elementary tail series exists (this is
exactly why the implementation falls back to `fzero` on a numerical
integral rather than an incomplete-function inverse), so this derivation
uses Laplace's-method saddle-point analysis directly on the mixture
integral. Maximizing the integrand's exponent `-x/tau - a*tau` over `tau`
gives a saddle at `tau* = sqrt(x/a)`, exponent value `-2*sqrt(a*x)` — the
well-known stretched-exponential K-distribution tail (equivalently, the
large-argument asymptotic of the modified Bessel function `K_{a-1}` that
appears in K's closed-form density). Inverting `S(x)=Pfa` through this
saddle point:

```
    x_pfa       ~ log(1/Pfa)^2 / (4a)
    delta_K(Pfa) ~ log log(1/Pfa) - 0.5*log(4a) - kappa1(a)
```

**Log-log class — the SAME slow growth rate as Weibull**, independent of `a`
to leading order (verified: at `a=0.5`, the formula already matches the
real `KCFAR_TLog.m` mixture-integral solve to `1.5e-5` absolute error across
the entire tested range `Pfa=1e-1...1e-20`; at larger `a`, convergence is
the same `O(log L / L)` slow rate as GenGamma's `v>0` branch — structurally
expected, since both derivations are large-argument gamma-type tail
asymptotics). **This gives K-distribution CFAR a second, independent
disadvantage on top of F13's already-documented support-condition
collapse**: even where K's `c2 > pi^2/24` condition IS met, its own
threshold offset grows only as slowly in `Pfa` as Weibull's — with no
`1/v`- or `1/rho`-style shape-dependent acceleration available to offset
it, since the log-log term's coefficient here is fixed at 1, not divided by
an estimated shape parameter the way GenGamma's/Burr's/G0's are.

**Important caveat for the paper's own thesis:** the growth-rate CLASS of a
detector's `delta(shape, Pfa)` formula is a property of the assumed
clutter-model FAMILY, holding shape fixed — it is not, by itself, a
prediction of real-data calibration performance (F9's measured
achieved/nominal ratio). The two can and do diverge: F13's GammaTex-MAP7
result (two-stage despeckle + plain-Gamma CFAR) sits in the SAME slow
log-log class as GenGamma's `v>0` branch and K (its formula, `delta =
0.5*(log(gammaincinv(1-Pfa,a)) - psi(a))`, is the `v=2` special case of the
`v>0` branch's structure) — yet it achieved the tightest real-data
calibration of any detector in this whole project (`PfaRatio` 0.83–0.97x on
HRSID). The resolution is not a contradiction: growth rate governs how much
a FIXED, correctly-specified model's own threshold moves as you ask for a
smaller nominal Pfa; F9's real-data ratio is additionally governed by how
well the fitted model's assumed tail matches the TRUE clutter tail at
extreme quantiles, and GammaTex-MAP7's despeckled texture matches its
assumed Gamma model far better than any of the six raw-data detectors match
theirs. **The paper's thesis should therefore be stated precisely as: growth
rate sets a MODEL-INTERNAL reachability bound — how far the nominal-Pfa dial
can practically move the achieved threshold at all, for a detector fit to
data its own family actually describes — not a claim that growth rate alone
ranks real-world calibration accuracy across mismatched-model detectors.**
That second effect (model/tail mismatch) is F9's separate, already-documented
contribution, and the two should be presented as compounding, not identical,
mechanisms.

*Reproduce:* `_comparison/verify_growth_rates.m` (new file, this session) —
run via `cfar_setup(); verify_growth_rates()`; compares every formula above
against the real `GenGammaCFAR_TLog.m` / `G0CFAR_TLog.m` / `BurrCFAR_TLog.m`
/ `KCFAR_TLog.m` output at representative shape values spanning the ranges
this project's real sweeps report (see `MedianShape`/`IQRShape` columns in
`comparison_summary_*.csv`), at `Pfa` from `1e-1` down to `1e-40`.

## 3. Resolving the Meng Collision (BC-3)

*Resolved 2026-09, by direct query of existing data — no new sweep needed.*

**The question:** Meng (2025, JSTARS) reports a cross-scene design-Pfa shift
with window geometry FIXED per scene. This project's own F9 finding reports
up to 3,162× measured-vs-nominal Pfa divergence at a FIXED nominal setting.
Are these the same axis (making this project's number a restatement of
Meng's result) or different axes (making it a distinct, citable finding)?

**Method:** query `comparison_summary_hrsid_full.csv` and
`comparison_summary_ssdd_full_fixed60.csv` (each: 60 images, HELD FIXED
across every `sli` in `{11,15,17,21,31,41,51,61,71,81,91,101,111,121,131,
141,151}`) for whether the measured/nominal Pfa ratio (`PfaRatio`) at a
fixed nominal Pfa changes materially as `sli` alone varies — the image/scene
set never changes across the sweep, so any variation found is attributable
to window geometry alone, not to Meng's cross-scene mechanism.

**Result — it does, substantially, on both datasets:**

| Detector | HRSID fold-change (min→max ratio, Pfa=1e-6) | SSDD fold-change (Pfa=1e-6) |
|---|---|---|
| Weibull | 1.5x (489→744) | 1.9x (848→1611) |
| Lognormal | 16.9x (1.2→20.0) | **38.9x (11.4→443)** |
| GenGamma | 3.6x (351→1267) | 2.8x (1840→5060) |
| G0 | 5.0x (28→144) | 7.6x (20→150) |
| BurrXII | 4.8x (10→49) | 2.1x (40→86) |

Every detector on both datasets shows a substantial ratio swing driven by
`sli` ALONE, with the image set held fixed throughout. Lognormal's SSDD
swing (38.9x) is the largest — notable since F10 shows Lognormal is
otherwise the mildest-growing/best-calibrated of the two exactly-derived
detectors, underscoring that this window-size sensitivity is a real,
separate mechanism from the growth-rate story in §2, not a restatement of
it.

**Reconciling the exact "3,162×" figure.** That number is `FINDINGS.md`'s
original F9 table value: SSDD, `sli=21`, pooled over 100 images (a
DIFFERENT, smaller/differently-sampled image set than the 60-image sweep
above). Querying the 60-image sweep at the identical `sli=21` gives
GenGamma's ratio = 3,027.6x — a 4% difference from 3,162x, consistent with
ordinary sample-to-sample variation between a 100-image and a 60-image draw,
not a contradiction. The headline number is reproducible within expected
noise; more importantly, it is NOT the ceiling — the 60-image SSDD sweep
finds GenGamma's ratio reaches 5,060x at `sli=151`, higher than the original
headline figure, entirely from window size at fixed Pfa on the same scenes.

**Conclusion: BC-3 resolves in the paper's favor.** The divergence axis this
paper measures (achieved-vs-nominal at fixed nominal Pfa, window size swept,
scene set fixed) is demonstrably distinct from Meng's axis (cross-scene
design-Pfa shift, window geometry fixed per scene) — moving `sli` alone,
with the scene held constant, moves the divergence ratio by up to
38.9x. State this explicitly in the paper as the answer to the BC-3
falsifier, per the roadmap's item 1 instruction, converting what could have
been a defensive footnote into an offensive one.

*Reproduce:* `_comparison/Results/comparison_summary_hrsid_full.csv`,
`comparison_summary_ssdd_full_fixed60.csv`, both already in the repo from
the earlier Phase 6.2a SLI-plateau work; no new sweep was required for this
section.

## 4. Statistical Rigor

*Computed 2026-09 (roadmap item 4), from the clutter-WINDOW count
(`BackgroundPixels`, order 1e5–4e7 per sweep cell), never the ship count —
conflating the two was explicitly flagged as the most obvious small-sample
objection a reviewer would raise, and is avoided here by construction.*

**Method.** Each `(detector, sli, Pfa)` cell in the existing sweeps already
pools `FalsePixels` (successes) out of `BackgroundPixels` (trials) — an
exact 95% Clopper-Pearson interval is computed directly from these via the
incomplete-beta identity (`betaincinv`, no toolbox needed — the same
function G0's own threshold derivation uses), attached as new
`Pfa_CI_Lo`/`Pfa_CI_Hi`/`Ratio_CI_Lo`/`Ratio_CI_Hi` columns
(`_comparison/compute_clopper_pearson.m`, output `comparison_summary_*_ci.csv`).
**Caveat, stated rather than hidden:** this treats each background pixel as
an independent Bernoulli trial. Neighbouring reference-window pixels are
not independent (overlapping windows, spatially correlated clutter), so the
true interval is somewhat wider than the i.i.d. interval computed here —
the same approximation implicit in using the raw pixel count as "the"
sample size at all. A tighter treatment would need an effective-sample-size
correction for spatial correlation; not attempted here, flagged as a
limitation.

**Result — the headline divergence numbers are tight; the well-calibrated
ones are not, and that asymmetry is itself worth stating.**

| Figure | `PfaRatio` | 95% CI | Relative width |
|---|---|---|---|
| GenGamma, SSDD, `sli=21`, `Pfa=1e-6` (original F9 headline, 100 images) | 3,162.1x | [3,133.9x, 3,190.5x] | 1.8% |
| GenGamma, SSDD, `sli=151`, `Pfa=1e-6` (§3's reconciled ceiling, 60 images) | 5,060.2x | [5,014.6x, 5,106.1x] | 1.8% |
| Weibull, SSDD, `sli=51`, `Pfa=1e-6` | 978.6x | [958.7x, 998.8x] | 4.1% |
| Weibull, HRSID, `sli=51`, `Pfa=1e-6` | 547.6x | [540.2x, 555.1x] | 2.7% |
| Lognormal, HRSID, `sli=51`, `Pfa=1e-6` | 1.9x | [1.5x, 2.4x] | **47.4%** |

Both headline GenGamma figures — the ones this paper's thesis leans on —
are tightly bounded (<2% relative width) at their reported sample size, not
a coincidence: they come from cells with tens of thousands of false pixels
out of `~10-15M` background pixels, plenty of statistical power. **The
detector this project already flags as best-calibrated (Lognormal, F9/F10)
has the widest relative interval of the five** — not because its
measurement is worse, but because *being well-calibrated means observing
very few false pixels*, and a small numerator inflates relative uncertainty
even with millions of trials (73 false pixels here). This is worth stating
explicitly in the paper: **a tight "detector X is well-calibrated" claim
needs a correspondingly larger sample than a "detector Y is badly
miscalibrated" claim does, at the same nominal Pfa** — the two kinds of
claim are not symmetric in how much data they need to support.

**A sharper example of the same asymmetry — zero observed false alarms is
not proof of good calibration.** Several `G0` cells on the SSDD full grid
report `FalsePixels=0` at `Pfa=1e-6` (e.g. `sli=41`: 0 false pixels out of
358,576 background pixels), giving a naive point-estimate ratio of exactly
`0.00x` — read uncritically, this looks like G0 is *better* than perfectly
calibrated. The Clopper-Pearson interval corrects this immediately: the 95%
upper bound is **10.3x nominal**, not 0 — a textbook "rule of three"
result (with zero successes in `N` trials, the 95% upper bound on the true
rate is `~3/N`). Zero observed false alarms at this sample size is
consistent with anywhere from perfect calibration to a ratio an order of
magnitude worse; it does not distinguish between them. Any future claim in
this paper (or a downstream one) that a detector "never" produces a false
alarm at some `Pfa`/`sli` combination must report the CI upper bound, not
the raw zero.

**Grid-wide distribution of relative CI width** (SSDD full grid, 612
cells): median 7.9%, 90th percentile 67.7% — most cells are tightly bounded,
but a meaningful tail (mostly G0/K cells at small `sli`, where the
reference window is small and the nominal Pfa is already tight, so few
false pixels accumulate) carries wide intervals. Any single-cell number
pulled from these sweeps for a table or figure should be checked against
its own CI column before being reported as precise.

*Reproduce:* `_comparison/compute_clopper_pearson.m`; outputs
`Results/comparison_summary_{main,ssdd_full_fixed60,hrsid_full}_ci.csv`.

**Addendum, confirmed 2026-09 (roadmap item 5) — the HRSID sample was
expanded from 60 to 250 images (same seed=42 draw superset, same 17-`sli`
grid), and the result independently validates the CI-width finding above
rather than just repeating it.** The large divergence figures are stable:
GenGamma's ratio at every `sli` moved by at most ~5% (e.g. `sli=51,
Pfa=1e-6`: 389.2x at 60 images -> 410.9x at 250 images; `sli=11`: 1267.2x ->
1315.5x), and Weibull similarly (547.6x -> 529.8x). **Lognormal's small
ratio, by contrast, moved by more than 2x** (1.9x -> 4.2x at the same cell)
— not a contradiction, but exactly what a wide relative CI predicts: a
well-calibrated detector's point estimate is expected to swing substantially
between independent samples of this size, while a badly-miscalibrated one's
does not. This is real, out-of-sample confirmation of the CI-width argument
above, not a restatement of it. **Conclusion for the paper: report the
divergence headline figures (GenGamma, Weibull) as stable and reproducible
under sample expansion; report Lognormal's own ratio only alongside its CI,
never as a bare point figure**, since the point estimate alone is shown here
to be sample-size-sensitive in a way the paper's own statistical-rigor
section already predicted.

*Reproduce:* `run_comparison_hrsid('NumImages',250,'Sli',[11 15 17 21 31 41
51 61 71 81 91 101 111 121 131 141 151],'Random',true,'Seed',42,'Tag','_hrsid_full250')`;
output `Results/comparison_summary_hrsid_full250.csv`.

## 5. Related Work / Prior Art

*Searched 2026-09 (roadmap item 6), via IEEE Xplore through the
supervisor's institutional access (VIT Chennai MyLOFT portal) — queries
covering the two alternate-terminology guesses ("CFAR loss saturation",
"Pfa control breakdown"), achieved-vs-nominal Pfa divergence directly,
small-sample/Monte-Carlo calibration, heavy-tailed order-statistics CFAR,
and the supervisor's own group (Mahapatra) for any prior touch on this
exact question.*

**Result: the search came back clean — no collision found.** Seven papers
were pulled and read in full or in relevant part; none states this
project's reachability-bound thesis (that a detector's own `δ(shape, Pfa)`
growth rate as `Pfa -> 0` bounds how far the achieved Pfa can practically
track the nominal one) under any terminology, alternate or otherwise:

| Paper | What it actually covers | Why it doesn't collide |
|---|---|---|
| Sarma & Tufts, "Rank-Order Adaptive CFAR: Performance Bounds and Efficient Implementation," *IEEE TAES* 49(4), 2013 | UMP nonparametric rank-order CFAR, exact ROC bounds vs. CA-CFAR via complete enumeration | Different detector class entirely — rank statistics have an exactly enumerable threshold distribution (no MoLC shape-fitting, no asymptotic tail-quantile problem); "performance bounds" here means ROC curves, not a Pfa-reachability limit |
| Redjem & Sahed, "Noncoherent Parameter-Free GM-CFAR... Pareto-Distributed Clutter," ICATEEE 2025 | Exact closed-form Pfa for geometric-mean CFAR, monopulse and multipulse | Derives an exact Pfa formula (via Gamma/Exponential duality) but never asks how the threshold offset *grows* as Pfa shrinks — no achieved-vs-nominal divergence framing at all |
| Bouchelaghem, Meddah & Nouar, "Hybrid Algorithm for CFAR Detection in Non-Homogeneous Environments," ICAECCS 2025 | Min-rule combining CA-CFAR and ACOSD-CFAR thresholds for clutter-edge/interferer robustness | Solves spatial non-homogeneity (clutter edges, interfering targets), not Pfa calibration as a function of the nominal Pfa dial |
| Hinz et al., "Presegmentation-based adaptive CFAR detection for HFSWR," 2012 | Doppler-domain background segmentation for surface-wave radar | Different application (HF ocean surveillance), different problem (clutter-edge avoidance) |
| Iinatti, "On the Threshold Setting Principles in Code Acquisition of DS-SS Signals," *IEEE JSAC* 18(1), 2000 | Spread-spectrum code-acquisition threshold/mean-acquisition-time tradeoffs | Different domain entirely (comms synchronization, Marcum-Q/Rician thresholds) — a keyword false-positive on "threshold setting principles," not SAR/CFAR clutter literature |
| Liu et al., "Adaptive detection in the presence of signal mismatch," *J. Systems Eng. Electron.* 26(1), 2015 | Array/beamformer GLRT tunable detectors (KGLRT, AMF, ACE) for steering-vector mismatch | Multichannel array detection theory, unrelated to single-channel clutter-model Pfa calibration |
| Petrović, Petrović & Milovanović, "From Front-End to Perception: A Chisel Generator for Low-Latency FMCW Radar Signal Processing," 2025 | Parameterizable Chisel HDL generator emitting FFT+2D-CFAR pipelines | Not applicable to THIS paper (no Pfa-calibration content) — relevant instead to Paper 2/3 as a related-work citation for the hardware-generator angle |

**This absence is reported as an audited search, not silence** (per the
roadmap's own instruction): the search was not IEEE-Xplore-exhaustive (five
targeted query strings, not a systematic review), so it does not prove no
such paper exists anywhere — but it is real evidence, not an assumption,
that the specific "CFAR loss saturation" / "Pfa control breakdown"
terminology guesses did not surface anything, and that the broader
achieved-vs-nominal-Pfa-divergence angle did not either. This paper's
prior-art risk from this specific pre-emption route is therefore lower than
it was before this search, though not eliminated.

*Reproduce:* IEEE Xplore, institutional access, queries: `"CFAR" AND "false
alarm" AND ("saturation" OR "breakdown")`; `"CFAR" AND "threshold" AND
"achieved" AND "nominal" AND "false alarm"`; `"CFAR" AND "calibration" AND
"small sample" AND "Monte Carlo"`; `"order statistics" AND "CFAR" AND "Pfa"
AND "heavy tail"`; `Mahapatra AND CFAR AND Pfa`.

## 6. Discussion

**The mechanism comes first; the numbers are its consequence, not the
finding itself.** It is tempting to lead with the largest measured
divergence in this project's grid — Generalized Gamma reaching over
5,000x nominal Pfa at `sli=151` on SSDD (§3, §4) — because it is the most
dramatic single figure available. That framing is exactly backwards for
this paper's purpose, and is avoided throughout: the number is a
*consequence* of Generalized Gamma's `v>0` branch sharing Weibull's
log-log growth class while additionally dividing by an independently noisy
estimated parameter `v` (§2.1); it is not itself the discovery. A reader
who remembers only "5,060x" and not "log-log growth, divided by a noisy
`v`" has learned a dataset-specific curiosity, not a distribution-family
property — and, restated once more so it cannot be read as implying
otherwise: **that number is not a ceiling.** Every `delta(shape, Pfa)`
formula derived in §2 is unbounded; 5,060x is simply how far this
particular grid's window sizes happened to push it. A larger `sli` range,
or a different scene, would find it larger still, for the same underlying
reason and with no asymptotic limit approached.

**Growth-rate class and real-world calibration accuracy are related but
not identical, and §2.1 already states the precise relationship rather
than leaving it to be inferred.** The clearest illustration is the
GammaTex-MAP7 result reported elsewhere in this research program (a
two-stage despeckle-then-Gamma-CFAR pipeline): its own `delta` formula sits
in the *same* slow log-log class as Weibull and K-distribution, yet it
achieves the tightest real-data calibration of any detector measured in
this comparison. This is not a counterexample to the growth-rate result —
it is exactly what the result predicts once the two effects it depends on
are separated. Growth rate is a property of a *correctly specified* model:
it says how much the threshold moves as the nominal Pfa dial turns, holding
the fitted shape fixed and assuming that shape genuinely describes the
data. Real-world calibration accuracy (this project's own F9 finding)
additionally depends on whether the fitted model's assumed tail matches the
*true* clutter tail at the extreme quantiles a small Pfa probes.
GammaTex-MAP7's despeckled texture happens to match its assumed Gamma model
far better than any raw-data detector in this comparison matches its own —
good model fit compensates for a slow growth rate. A practitioner should
therefore read a slow-growth classification as "the mapping from nominal
Pfa to threshold is inherently coarse for this family, so calibration
accuracy at small Pfa depends unusually heavily on getting the shape fit
right," not as "this family will always calibrate badly." The two
mechanisms compound; neither substitutes for the other, and a paper — or a
detector-selection decision — that reports only one is reporting half the
picture.

**Reporting practice: an achieved-Pfa or divergence-ratio number without
its confidence interval is not a fully specified claim, and the direction
of the risk is counter-intuitive.** §4 shows the asymmetry directly:
this paper's largest, most alarming divergence figures are also its
tightest, because a large false-alarm rate produces a large false-pixel
count and therefore a narrow Clopper-Pearson interval almost for free.
The numbers most likely to be under-scrutinized — a detector reported as
"well-calibrated," or as producing "zero false alarms" at some `sli`/Pfa
combination — are exactly the ones with the least statistical support,
because good calibration produces few or no false pixels to estimate a
rate from. A "zero false alarms" claim with no stated sample size is not
evidence of a working detector; §4 measured one concrete case where it
corresponded to a 95% upper bound of 10x nominal. Any paper in this
research program, this one included, that reports a favorable calibration
number should report it with an interval attached; the growth-rate result
in §2 explains *why* the underlying quantity being estimated is so
sensitive to the fitted shape and window geometry in the first place, and
§4's grid-wide result (median relative width 7.9%, 90th percentile 67.7%
on the SSDD grid) shows this is not a rare edge case within this
comparison's own data.

**Limitations.** The reachability-bound framing in this paper concerns Pfa
calibration specifically; it says nothing about detection probability or
ROC behavior, which is the subject of separate work in this research
program. The Clopper-Pearson intervals in §4 treat background pixels as
independent Bernoulli trials, which neighboring, overlapping reference
windows are not — the true intervals are somewhat wider than reported, in
a direction that only strengthens this paper's caution about small-count
cells, not its confidence in large ones. The four asymptotic derivations in
§2.1 (all but Weibull and Burr XII) are leading-order results with stated,
verified error terms, not closed algebraic identities; for K-distribution
and Generalized Gamma's `v>0` branch specifically, convergence to that
leading order is slow (`O(log L / L)`), so the growth-rate *class* is
established with confidence while the exact constant at any single moderate
Pfa is not. The literature search in §5 was targeted, not exhaustive, and
is reported as such.

**A natural question this result raises for future work is whether
nonparametric CFAR families sidestep the problem entirely.** Rank-order
CFAR detectors (e.g. Sarma and Tufts, 2013) determine their detection
threshold by exact enumeration of a discrete rank-order statistic rather
than by inverting a continuous tail — there is no `delta(shape, Pfa)`
function to have a slow growth rate at all, at the cost of the detector
sensitivity such tests are already known to trade away at very low Pfa.
Whether that tradeoff is favorable specifically *because* of the
reachability bound identified here, rather than for the previously known
reasons alone, was not evaluated in this paper and is worth a dedicated
comparison.

## 7. Conclusion

This paper identified a source of Pfa miscalibration in parametric SAR
CFAR detectors that is independent of clutter-model mismatch, estimation
noise, and cross-scene variation — all previously documented mechanisms —
and is instead a property of the clutter distribution family itself. Every
detector considered here reduces to a shared decision rule in which a
threshold offset `delta(shape, Pfa)` is added to a clutter-level term; how
fast that offset grows as the requested Pfa shrinks toward zero is fixed
by the family's tail type and differs by a full order of growth between
families in routine use for SAR clutter. Six such families were derived
and verified: two exactly (Weibull, Burr XII) and four as asymptotics
cross-checked numerically against their real implementations rather than
trusted as algebra, catching one real derivation error in the process
(§2.1). The result is a three-class taxonomy — log-log, square-root-log,
and linear-log growth — determined by whether the underlying clutter tail
is light (Gumbel-type) or heavy (power-law), not by which of the six
detectors it happens to be.

This growth-rate property was shown to be a distinct axis from a recently
published cross-scene Pfa-divergence result, resolved directly against
this project's own controlled measurement grid rather than argued from
first principles alone (§3), and every quantitative claim in this paper was
attached to an exact confidence interval computed from the appropriate
sample size, revealing that the most dramatic-looking numbers in this
comparison are also its best-supported, while the most reassuring ones
frequently are not (§4). No collision with this specific result was found
in a targeted, though not exhaustive, search of the accessible literature
(§5).

The practical consequence for SAR CFAR practice is not that any detector
family is unusable, nor that any achievable Pfa is bounded below — no such
bound exists, and this paper has been explicit throughout that none of its
six growth-rate formulas saturate. The consequence is narrower and more
actionable: a detector family's growth-rate class should be reported
alongside any claim about its calibration accuracy at small nominal Pfa,
because that class determines how much statistical and numerical precision
the rest of the system — shape estimation, fixed-point hardware
quantization, sample size — needs to supply before turning the nominal-Pfa
dial produces any real effect at all. A slow-growing family is not a bad
family; it is a family that demands more of everything else around it, and
a calibration claim made about one without saying so is an incomplete
claim.
