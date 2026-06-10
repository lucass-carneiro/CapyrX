# CapyrX_WaveToy Code Review

**Reviewed:** 2026-06-10  
**Thorn:** `arrangements/CapyrX/CapyrX_WaveToy`  
**Scope:** All source files and parameter files in the thorn.

---

## Summary

The thorn implements a first-order reduction of the scalar wave equation in curved spacetime,
following Eqs. (56)–(58) of [gr-qc/0507004](https://arxiv.org/abs/gr-qc/0507004). The RHS is
written in flux-conservative form using global (physical-coordinate) derivatives obtained by
contracting local patch derivatives with the multipatch Jacobian.

Five bugs ranging from simulation-breaking to minor were found, plus several par-file defects.

---

## Bug 1 — CRITICAL: Missing `z²` term in `quad_gaussian::phi`

**File:** `src/quad_gaussian.hxx`, line 16–21

**Code:**
```cpp
return (A * (x - x0 + y - y0) * (x - x0 - y + y0) *
        (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
       (exp(...) *
        (pow(x - x0, 2) + pow(y - y0, 2)) *
        (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)));
```

The numerator factor intended to be `(z - z0)²` is written as `-2·z·z0 + z0²`, which equals
`(z - z0)² - z²`. The missing `z²` makes the factor equal to `(x-x0)² + (y-y0)²` only when
`z0 = 0`, which happens to cancel the matching denominator factor and produces the correct
quadrupolar mode `∝ sin²θ cos(2φ) / r²` in that special case. For **any nonzero `z0`** the
formula is wrong.

The same incorrect factor `- 2*z*z0 + pow(z0, 2)` propagates into all derivative functions
(`Dx`, `Dy`, `Dz`, `dDxdx`, `dDydy`, `dDzdz`), which are self-consistent with the buggy `phi`
but also incorrect.

**Fix:** Replace every occurrence of `-2 * z * z0 + pow(z0, 2)` with `pow(z - z0, 2)` in
`quad_gaussian.hxx`.

---

## Bug 2 — CRITICAL: Broken parameter assignments in par files

**Files:** `par/multipole.par:72–73`, `par/overlap.par:68`, `par/quad_gauss_cubed_sphere.par:60`

Three par files contain a bare `CapyrX_WaveToy:: = <value>` assignment — the parameter name is
missing entirely:

```
# multipole.par:72
CapyrX_WaveToy::    = $r0 / 2.5          # parameter name missing
CapyrX_WaveToy::dissipation_epsilon = $diss_eps  # parameter does not exist
```

```
# overlap.par:68
CapyrX_WaveToy::  = $r0 / 2.0            # parameter name missing

# quad_gauss_cubed_sphere.par:60
CapyrX_WaveToy::    = $r0 / 4.5          # parameter name missing
```

From context these lines should set `gaussian_R0`. The Cactus parameter parser will either abort
or silently ignore them, meaning the Gaussian radius is never set (remains at the default `1.0`).

Additionally, `multipole.par:73` uses `CapyrX_WaveToy::dissipation_epsilon` which does not exist
in `param.ccl`. The correct parameters are `cart_diss_eps` and `curv_diss_eps`.

**Fix:**
```
# multipole.par
CapyrX_WaveToy::gaussian_R0 = $r0 / 2.5
CapyrX_WaveToy::cart_diss_eps = $diss_eps
CapyrX_WaveToy::curv_diss_eps = $diss_eps

# overlap.par
CapyrX_WaveToy::gaussian_R0 = $r0 / 2.0

# quad_gauss_cubed_sphere.par
CapyrX_WaveToy::gaussian_R0 = $r0 / 4.5
```

---

## Bug 3 — CRITICAL: `sw_cubed_sphere.par` omits `ADMBaseX`

**File:** `par/sw_cubed_sphere.par`, line 1–9

```
ActiveThorns = "
    CarpetX IOUtil TimerReport ODESolvers
    CoordinatesX CapyrX_MultiPatch CapyrX_WaveToy
"
```

`ADMBaseX` is absent. The RHS reads `ADMBaseX::metric`, `lapse`, and `shift`, and the schedule
declares `READS: ADMBaseX::metric(everywhere)`. Without the thorn active, the metric grid
functions are uninitialized. Because the par file also sets `CarpetX::poison_undefined_values =
yes`, `gxx = gyy = gzz = NaN`, so `calc_det(g_dd)` returns NaN immediately and every evolved
variable becomes NaN at the first time step.

**Fix:** Add `ADMBaseX` to `ActiveThorns` and add a flat-metric initializer (e.g.,
`ADMBaseX::initial_data_setup_method = "init_some_levels"` with a suitable thorn, or use
`TmunuBaseX`'s flat initializer).

---

## Bug 4 — MODERATE: `phi_rhs` and `D_i_rhs` drop shift advection

**File:** `src/rhs.cxx`, lines 165, 172–174

```cpp
phi_rhs(p.I) = alpha_Pi(0, p.I);          // = α·Π only

Dx_rhs(p.I) = g_dalphaPi.dx;              // = ∂_x(α·Π) only
Dy_rhs(p.I) = g_dalphaPi.dy;
Dz_rhs(p.I) = g_dalphaPi.dz;
```

The correct evolution equations are:
```
∂_t φ = α Π + βⁱ Dᵢ
∂_t Dᵢ = ∂ᵢ(α Π + βʲ Dⱼ)
```

The shift advection term `βⁱ Dᵢ` is absent from both. For a flat Minkowski background
(`β = 0`) this is harmless. For Kerr-Schild backgrounds (used in `schw_wave.par`) the shift is
`O(M/r)` and these terms are not negligible. The code comment documents the limitation ("not
valid if the metric is not static") but the `schw_wave.par` parameter file silently enables this
thorn with a non-trivial shift, giving quietly wrong results for the φ and D_i time series.

The comment also notes that the extrinsic curvature source term `-α K Π` in the Π equation is
absent; this is a separate but related limitation.

**Consequence:** For the Kerr-Schild run the φ field drifts on the orbital/light-crossing
timescale compared to the correct evolution. The energy norm `||Π||² + g^{ij} D_i D_j` will not
be conserved even approximately.

---

## Bug 5 — MODERATE: `error.cxx` aborts for Gaussian initial conditions

**File:** `src/error.cxx`, lines 53–55

```cpp
} else {
  CCTK_ERROR("Unknown initial condition");
}
```

The error function is only scheduled under the condition `CCTK_EQUALS(initial_condition,
"standing wave")` (schedule.ccl line 75), so the `else` branch cannot be reached through the
normal schedule. However, any future change relaxing the schedule condition, or a direct call
during debugging, will cause an immediate abort for the "Gaussian" and "Quadrupolar Gaussian"
initial conditions, which do have well-defined analytical solutions or at least energy norms that
could serve as error measures.

**Fix:** Replace `CCTK_ERROR` with a no-op or a per-condition error computation.

---

## Bug 6 — MINOR: `NewRadX_Apply` uses the same `rad_power` for all variables

**File:** `src/rhs.cxx`, lines 230–243

```cpp
NewRadX_Apply(cctkGH, phi,  phi_rhs,  ..., 0.0, 1.0, rad_power);
NewRadX_Apply(cctkGH, Pi,   Pi_rhs,   ..., 0.0, 1.0, rad_power);
NewRadX_Apply(cctkGH, Dx,   Dx_rhs,   ..., 0.0, 1.0, rad_power);
...
```

The second-order sibling `CapyrX_WaveToy_SecondOrder` (rhs.cxx line 71) correctly uses
`rad_power + 1.0` for `rho_rhs` because the conjugate momentum falls off one power faster than
the field. Here, φ ∼ r⁻¹ and Π, Dᵢ ∼ r⁻², so `Pi_rhs`, `Dx_rhs`, `Dy_rhs`, `Dz_rhs` should
use `rad_power + 1.0`. Using the same power for all variables causes the radiative boundary
condition to apply the wrong falloff to the momentum variables, degrading wave extraction
accuracy at finite radii.

---

## Bug 7 — MINOR: `quad_gauss::Pi` returns `int 0` instead of `T{0}`

**File:** `src/quad_gaussian.hxx`, line 28

```cpp
return 0;   // should be T{0}
```

For a template function returning type `T`, returning the integer literal `0` relies on implicit
conversion. If `T` is ever instantiated with a type that is not implicitly constructible from
`int` (e.g., an interval-arithmetic or dual-number type used for code verification), this will
fail to compile. The fix is `return T{0};`.

---

## Bug 8 — MINOR: Dead code in `standing_wave.hxx` and `quad_gaussian.hxx`

**Files:** `src/standing_wave.hxx` lines 62–120, `src/quad_gaussian.hxx` lines 171–379

Both headers export `dPidx`, `dPidy`, `dPidz`, `dDxdx`, `dDydy`, `dDzdz` — spatial derivatives
of the analytical solution. These functions are never called anywhere in the thorn (neither in
`error.cxx` nor in any test harness). The dead code adds maintenance burden: the
`quad_gaussian.hxx` versions are hundreds of lines long and embed the same `z0` bug as `phi`.

---

## Bug 9 — MINOR: `use_newradx` scheduled AFTER `CapyrX_WaveToy_Dissipation`

**File:** `schedule.ccl`, lines 51–59

```
SCHEDULE CapyrX_WaveToy_ApplyNewRadX IN ODESolvers_RHS
  AFTER CapyrX_WaveToy_Dissipation
```

If `add_dissipation = no`, `CapyrX_WaveToy_Dissipation` is not registered, so the `AFTER`
ordering constraint references a non-existent item. Cactus treats a missing `AFTER` target as a
no-op, so this is not a runtime crash, but it documents a false dependency. The NewRadX routine
should be scheduled `AFTER CapyrX_WaveToy_RHS` (and optionally also after
`CapyrX_WaveToy_Dissipation` with a conditional guard), or the two conditionals should be nested.

---

## Minor Code Quality Issues

| Location | Issue |
|---|---|
| `rhs.cxx:19` | Variable named `den` holds `1/(12·DX)` — i.e. the *reciprocal* of the denominator. Misleading naming. |
| `energy.cxx:18` | Lambda marked `CCTK_HOST CCTK_DEVICE` while all other device lambdas in the thorn use `CCTK_DEVICE` only. Inconsistent; could affect GPU dispatch. |
| `rhs.cxx` vs `local_derivatives.hxx` | The `diss_5` stencil is duplicated verbatim in both thorns. Should be shared via a header. |
| `rhs.cxx` | Three separate functions `c4o_1_0_0`, `c4o_0_1_0`, `c4o_0_0_1` could be a single `template<int dir>` as in `CapyrX_WaveToy_SecondOrder/src/local_derivatives.hxx`. |
| `interface.ccl:3` | Comment says implementation follows gr-qc/0512001v1 but `rhs.cxx` references gr-qc/0507004. The governing equations reference should be consistent. |
| `interface.ccl:21` | The `state` group tag `dependents="energy error"` lists `error` as a dependent of `state`. But the `error` group stores φ_error etc., not physical state. If the scheduler uses this tag for automatic recomputation triggering, listing `error` here may cause spurious recomputes. |
