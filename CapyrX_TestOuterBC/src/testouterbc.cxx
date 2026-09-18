#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>

#include <array>
#include <cmath>
#include <limits>

namespace CapyrX::TestOuterBC {

// Must live at namespace scope (not inside a function): nvcc forbids a
// function-local type from being used as the type of a variable captured by
// an extended __device__ lambda.
enum class SmoothKind { z_global, parabola, one_over_r };

extern "C" void CapyrX_TestOuterBC_write_nan_test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_TestOuterBC_write_nan_test;
  DECLARE_CCTK_PARAMETERS;

  const auto nan_value{std::numeric_limits<CCTK_REAL>::quiet_NaN()};

  const bool inject_ghost{CCTK_Equals(nan_injection_mode, "ghost") ||
                          CCTK_Equals(nan_injection_mode, "ghost_and_interior")};
  const bool inject_interior{
      CCTK_Equals(nan_injection_mode, "interior") ||
      CCTK_Equals(nan_injection_mode, "ghost_and_interior")};

  // Per-face flag: true where this patch's face is the genuine physical
  // outer boundary, i.e. not shared with another patch. Same call
  // write_color already makes to size its overlap band.
  std::array<CCTK_INT, 2 * Loop::dim> is_interpatch_face{};
  if (CCTK_IsFunctionAliased("MultiPatch_GetBoundarySpecification2"))
    MultiPatch_GetBoundarySpecification2(grid.patch, 2 * Loop::dim,
                                         is_interpatch_face.data());

  Loop::vect<bool, Loop::dim> is_outer_lo, is_outer_hi;
  for (int d = 0; d < Loop::dim; ++d) {
    is_outer_lo[d] = !is_interpatch_face[2 * d + 0];
    is_outer_hi[d] = !is_interpatch_face[2 * d + 1];
  }

  const auto gsh{grid.gsh};
  const auto lbnd{grid.lbnd};
  const auto nghostzones{grid.nghostzones};

  // Baseline fill and NaN injection in a single pass: each cell's final value
  // is fully determined by its own position (either the 1.0 sentinel or
  // NaN), so there is no need to fill first and overwrite afterwards.
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        bool inject = false;
        for (int d = 0; d < Loop::dim; ++d) {
          const auto gI{p.I[d] + lbnd[d]};
          if (is_outer_lo[d]) {
            if (inject_ghost && gI < nghostzones[d])
              inject = true;
            if (inject_interior && gI == nghostzones[d])
              inject = true;
          }
          if (is_outer_hi[d]) {
            if (inject_ghost && gI >= gsh[d] - nghostzones[d])
              inject = true;
            if (inject_interior && gI == gsh[d] - nghostzones[d] - 1)
              inject = true;
          }
        }
        nan_test(p.I) = inject ? nan_value : 1.0;
      });

  // Snapshot into nan_test_pre before returning, mirroring color_pre: "SYNC:
  // nan_test" runs immediately after this routine and overwrites nan_test's
  // ghost zones via interpatch interpolation, so nan_test_pre (never synced)
  // preserves exactly what was injected.
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        nan_test_pre(p.I) = nan_test(p.I);
      });
}

extern "C" void CapyrX_TestOuterBC_write_smooth_test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_TestOuterBC_write_smooth_test;
  DECLARE_CCTK_PARAMETERS;

  // Which analytic field to use, resolved once on the host so the device
  // lambda only carries an int (CCTK_Equals on the raw keyword string inside a
  // device kernel is neither needed nor cheap).
  SmoothKind kind;
  if (CCTK_Equals(smooth_field, "z_global"))
    kind = SmoothKind::z_global;
  else if (CCTK_Equals(smooth_field, "parabola"))
    kind = SmoothKind::parabola;
  else if (CCTK_Equals(smooth_field, "one_over_r"))
    kind = SmoothKind::one_over_r;
  else
    CCTK_VERROR("Unknown smooth_field type %s", smooth_field);

  // Fill EVERYWHERE (loop_all), including the outer-BC ghost cells: this is
  // what lets the boundary_* = none run keep exact analytic values in the
  // outer ghosts (measurement (c)). The vertex_coords field is valid in the
  // ghost zones too (write_color/write_nan_test already read it there), so the
  // exact value in an outer ghost uses that ghost's genuine global coordinate.
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        using std::sqrt;
        const auto x{vcoordx(p.I)};
        const auto y{vcoordy(p.I)};
        const auto z{vcoordz(p.I)};

        CCTK_REAL val;
        switch (kind) {
        case SmoothKind::z_global:
          val = z;
          break;
        case SmoothKind::parabola:
          val = x * x + y * y + z * z;
          break;
        case SmoothKind::one_over_r:
          val = 1.0 / sqrt(x * x + y * y + z * z);
          break;
        default:
          val = 0.0; // unreachable; keeps the compiler quiet about `val`
          break;
        }
        smooth(p.I) = val;
      });

  // Snapshot into smooth_pre before returning, mirroring nan_test_pre /
  // color_pre: "SYNC: smooth" runs immediately after this routine and
  // overwrites smooth's ghost zones (interpatch cells via interpolation,
  // outer-BC cells via the boundary kernel unless boundary_* = none), so
  // smooth_pre (never synced) preserves the exact analytic value everywhere.
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        smooth_pre(p.I) = smooth(p.I);
      });

  // BUGFIX_TODO.md C6 instrument (R2, gated off by default): overwrite this
  // patch's outer-BC ghost shell with smooth_sentinel, after the smooth_pre
  // snapshot above so the ground truth stays the exact analytic field. With
  // boundary_* = none the BC pass below is a no-op, so the sentinel survives
  // into SYNC and any neighbouring patch's interpatch interpolation that
  // reads this shell (a slaved fill included) carries it forward scaled by
  // its stencil weight -- see param.ccl.
  if (smooth_sentinel != 0.0) {
    std::array<CCTK_INT, 2 * Loop::dim> is_interpatch_face{};
    if (CCTK_IsFunctionAliased("MultiPatch_GetBoundarySpecification2"))
      MultiPatch_GetBoundarySpecification2(grid.patch, 2 * Loop::dim,
                                           is_interpatch_face.data());

    Loop::vect<bool, Loop::dim> is_outer_lo, is_outer_hi;
    for (int d = 0; d < Loop::dim; ++d) {
      is_outer_lo[d] = !is_interpatch_face[2 * d + 0];
      is_outer_hi[d] = !is_interpatch_face[2 * d + 1];
    }

    const auto gsh{grid.gsh};
    const auto lbnd{grid.lbnd};
    const auto nghostzones{grid.nghostzones};

    grid.loop_all_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
          bool is_outer_ghost = false;
          for (int d = 0; d < Loop::dim; ++d) {
            const auto gI{p.I[d] + lbnd[d]};
            if (is_outer_lo[d] && gI < nghostzones[d])
              is_outer_ghost = true;
            if (is_outer_hi[d] && gI >= gsh[d] - nghostzones[d])
              is_outer_ghost = true;
          }
          if (is_outer_ghost)
            smooth(p.I) = smooth_sentinel;
        });
  }
}

} // namespace CapyrX::TestOuterBC
