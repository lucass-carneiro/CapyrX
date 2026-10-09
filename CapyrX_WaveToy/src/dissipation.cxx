#include "cctk_Config.h"
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>

#include <AMReX_Arena.H>
#include <AMReX_GpuAtomic.H>

#include <cmath>
#include <limits>

namespace CapyrX::WaveToy {

using namespace Loop;

template <std::size_t dir>
static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
diss_5(const PointDesc &p,
       const GF3D2<const CCTK_REAL> &gf) noexcept -> CCTK_REAL {
  const auto fac{(1.0 / 64.0) * (1.0 / p.DX[dir])};
  const auto stencil{gf(p.I - 3 * p.DI[dir]) - 6.0 * gf(p.I - 2 * p.DI[dir]) +
                     15.0 * gf(p.I - 1 * p.DI[dir]) - 20.0 * gf(p.I) +
                     15.0 * gf(p.I + 1 * p.DI[dir]) -
                     6.0 * gf(p.I + 2 * p.DI[dir]) + gf(p.I + 3 * p.DI[dir])};
  return fac * stencil;
}

// diss_5<d> already carries 1/h_d, so weights[d] = R_d * h_d
static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
apply_ko_diss(const PointDesc &p, const CCTK_REAL (&weights)[3],
              const GF3D2<const CCTK_REAL> &gf) noexcept -> CCTK_REAL {
  return weights[0] * diss_5<0>(p, gf) + weights[1] * diss_5<1>(p, gf) +
         weights[2] * diss_5<2>(p, gf);
}

// Nyquist damping rates R_d of the KO operator on a curvilinear patch, and the
// local physical spacing h_min they are built from. h_min is NaN or <= 0 only
// for a degenerate coordinate map.
struct CurvKORates {
  CCTK_REAL rates[3];
  CCTK_REAL h_min;
};

CCTK_DEVICE CCTK_HOST static inline CurvKORates
curv_ko_rates(const PointDesc &p, const CCTK_REAL (&eps)[3],
              const GF3D2<const CCTK_REAL> &vJ_da_dx,
              const GF3D2<const CCTK_REAL> &vJ_da_dy,
              const GF3D2<const CCTK_REAL> &vJ_da_dz,
              const GF3D2<const CCTK_REAL> &vJ_db_dx,
              const GF3D2<const CCTK_REAL> &vJ_db_dy,
              const GF3D2<const CCTK_REAL> &vJ_db_dz,
              const GF3D2<const CCTK_REAL> &vJ_dc_dx,
              const GF3D2<const CCTK_REAL> &vJ_dc_dy,
              const GF3D2<const CCTK_REAL> &vJ_dc_dz) {
  using std::fmin, std::sqrt;

  const CCTK_REAL grad_norms[3] = {
      sqrt(vJ_da_dx(p.I) * vJ_da_dx(p.I) + vJ_da_dy(p.I) * vJ_da_dy(p.I) +
           vJ_da_dz(p.I) * vJ_da_dz(p.I)),

      sqrt(vJ_db_dx(p.I) * vJ_db_dx(p.I) + vJ_db_dy(p.I) * vJ_db_dy(p.I) +
           vJ_db_dz(p.I) * vJ_db_dz(p.I)),

      sqrt(vJ_dc_dx(p.I) * vJ_dc_dx(p.I) + vJ_dc_dy(p.I) * vJ_dc_dy(p.I) +
           vJ_dc_dz(p.I) * vJ_dc_dz(p.I))};

  const CCTK_REAL h_min =
      fmin(p.DX[0] / grad_norms[0],
           fmin(p.DX[1] / grad_norms[1], p.DX[2] / grad_norms[2]));

  return CurvKORates{{eps[0] / h_min, eps[1] / h_min, eps[2] / h_min}, h_min};
}

CCTK_DEVICE CCTK_HOST static inline bool valid_h_min(const CCTK_REAL h_min) {
  using std::isfinite;
  return isfinite(h_min) && h_min > 0.0;
}

// Factor (<= 1) that brings dt * sum_d R_d down to the cap
CCTK_DEVICE CCTK_HOST static inline CCTK_REAL
step_damping_scale(const CCTK_REAL step_damping, const CCTK_REAL cap) {
  return step_damping > cap ? cap / step_damping : 1.0;
}

// dt * sum_d R_d on a Cartesian patch, R_d = eps / h_d
static inline CCTK_REAL cartesian_step_damping(const GridDescBase &grid,
                                               const CCTK_REAL dt,
                                               const CCTK_REAL eps) {
  return dt * eps *
         (1.0 / grid.dx[0] + 1.0 / grid.dx[1] + 1.0 / grid.dx[2]);
}

static inline bool patch_is_cartesian(const int patch) {
  return static_cast<bool>(
             CCTK_IsFunctionAliased("MultiPatch_PatchIsCartesian"))
             ? static_cast<bool>(MultiPatch_PatchIsCartesian(patch))
             : true;
}

extern "C" void CapyrX_WaveToy_Dissipation(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_WaveToy_Dissipation;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const CCTK_REAL cap = diss_max_step_damping;

  if (patch_is_cartesian(grid.patch)) {
    // On a Cartesian patch R_d = cart_diss_eps / h_d, so the clamp is the same
    // at every point. Leaving cart_diss_eps untouched when it is not hit keeps
    // this path bitwise identical to plain KO.
    const CCTK_REAL eps =
        cart_diss_eps *
        step_damping_scale(cartesian_step_damping(grid, dt, cart_diss_eps),
                           cap);
    const CCTK_REAL weights[3] = {eps, eps, eps};

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          phi_rhs(p.I) += apply_ko_diss(p, weights, phi);
          Pi_rhs(p.I) += apply_ko_diss(p, weights, Pi);
          Dx_rhs(p.I) += apply_ko_diss(p, weights, Dx);
          Dy_rhs(p.I) += apply_ko_diss(p, weights, Dy);
          Dz_rhs(p.I) += apply_ko_diss(p, weights, Dz);
        });
    return;
  }

  const CCTK_REAL curv_eps[3] = {curv_diss_eps_a, curv_diss_eps_b,
                                 curv_diss_eps_c};

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const auto ko =
            curv_ko_rates(p, curv_eps, vJ_da_dx, vJ_da_dy, vJ_da_dz, vJ_db_dx,
                          vJ_db_dy, vJ_db_dz, vJ_dc_dx, vJ_dc_dy, vJ_dc_dz);

        // CapyrX_WaveToy_check_dissipation aborts on this at BASEGRID; this only
        // guards against a map that changed since.
        if (!valid_h_min(ko.h_min)) {
#if !defined(__CUDACC__) && !defined(__HIP_PLATFORM_AMD__) &&                   \
    !defined(__HIP_PLATFORM_HCC__) && !defined(__INTEL_LLVM_COMPILER)
          CCTK_VERROR("KO dissipation: degenerate coordinate map, local "
                      "physical spacing h_min = %g",
                      double(ko.h_min));
#endif
          amrex::Abort();
        }

        const CCTK_REAL scale = step_damping_scale(
            dt * (ko.rates[0] + ko.rates[1] + ko.rates[2]), cap);

        const CCTK_REAL weights[3] = {scale * ko.rates[0] * p.DX[0],
                                      scale * ko.rates[1] * p.DX[1],
                                      scale * ko.rates[2] * p.DX[2]};

        phi_rhs(p.I) += apply_ko_diss(p, weights, phi);
        Pi_rhs(p.I) += apply_ko_diss(p, weights, Pi);
        Dx_rhs(p.I) += apply_ko_diss(p, weights, Dx);
        Dy_rhs(p.I) += apply_ko_diss(p, weights, Dy);
        Dz_rhs(p.I) += apply_ko_diss(p, weights, Dz);
      });
}

namespace {
struct DissipationCheck {
  CCTK_REAL max_step_damping;      // largest dt * sum_d R_d at a clamped point
  CCTK_REAL min_h_min;             // smallest h_min at a clamped point
  CCTK_REAL nclamped;              // clamped points on curvilinear patches
  CCTK_REAL ninvalid;              // points with a degenerate coordinate map
  CCTK_REAL cart_max_step_damping; // largest clamped value, Cartesian patches
  bool dt_invalid;                 // dt was not usable, nothing was checked
};

DissipationCheck *dissipation_check_buffer() {
  static DissipationCheck *const buffer = static_cast<DissipationCheck *>(
      amrex::The_Managed_Arena()->alloc(sizeof(DissipationCheck)));
  return buffer;
}

// amrex::Gpu::Atomic::{Max,Min} are not atomic on the host, and CarpetX calls
// local routines from several OpenMP threads at once
CCTK_DEVICE CCTK_HOST inline void atomic_max(CCTK_REAL *const m,
                                             const CCTK_REAL value) {
#if AMREX_DEVICE_COMPILE
  amrex::Gpu::Atomic::Max(m, value);
#else
#pragma omp critical(CapyrX_WaveToy_dissipation_check)
  *m = *m > value ? *m : value;
#endif
}

CCTK_DEVICE CCTK_HOST inline void atomic_min(CCTK_REAL *const m,
                                             const CCTK_REAL value) {
#if AMREX_DEVICE_COMPILE
  amrex::Gpu::Atomic::Min(m, value);
#else
#pragma omp critical(CapyrX_WaveToy_dissipation_check)
  *m = *m < value ? *m : value;
#endif
}

} // namespace

extern "C" void CapyrX_WaveToy_reset_dissipation_check(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_WaveToy_reset_dissipation_check;

  auto *const buffer = dissipation_check_buffer();
  buffer->max_step_damping = 0.0;
  buffer->min_h_min = std::numeric_limits<CCTK_REAL>::infinity();
  buffer->nclamped = 0.0;
  buffer->ninvalid = 0.0;
  buffer->cart_max_step_damping = 0.0;
  buffer->dt_invalid = false;
}

extern "C" void CapyrX_WaveToy_check_dissipation(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_WaveToy_check_dissipation;
  DECLARE_CCTK_PARAMETERS;

  using std::isfinite;

  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const CCTK_REAL cap = diss_max_step_damping;
  const CCTK_REAL curv_eps[3] = {curv_diss_eps_a, curv_diss_eps_b,
                                 curv_diss_eps_c};
  auto *const buffer = dissipation_check_buffer();

  if (!(isfinite(dt) && dt > 0.0)) {
#pragma omp atomic write
    buffer->dt_invalid = true;
    return;
  }

  if (patch_is_cartesian(grid.patch)) {
    // Uniform over the patch: no loop needed
    const CCTK_REAL step_damping =
        cartesian_step_damping(grid, dt, cart_diss_eps);
    if (step_damping > cap)
      atomic_max(&buffer->cart_max_step_damping, step_damping);
    return;
  }

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const auto ko =
            curv_ko_rates(p, curv_eps, vJ_da_dx, vJ_da_dy, vJ_da_dz, vJ_db_dx,
                          vJ_db_dy, vJ_db_dz, vJ_dc_dx, vJ_dc_dy, vJ_dc_dz);

        if (!valid_h_min(ko.h_min)) {
          amrex::HostDevice::Atomic::Add(&buffer->ninvalid, CCTK_REAL(1));
          return;
        }

        const CCTK_REAL step_damping =
            dt * (ko.rates[0] + ko.rates[1] + ko.rates[2]);

        if (step_damping > cap) {
          amrex::HostDevice::Atomic::Add(&buffer->nclamped, CCTK_REAL(1));
          atomic_max(&buffer->max_step_damping, step_damping);
          atomic_min(&buffer->min_h_min, ko.h_min);
        }
      });
}

extern "C" void CapyrX_WaveToy_report_dissipation_check(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_WaveToy_report_dissipation_check;
  DECLARE_CCTK_PARAMETERS;

  // CarpetX synchronizes the device after every local routine, so the buffer is
  // complete here. Each process reports on its own points.
  const auto *const buffer = dissipation_check_buffer();

  if (buffer->dt_invalid) {
    CCTK_WARN(CCTK_WARN_ALERT,
              "KO dissipation: the time step is not set at this BASEGRID, so "
              "it was not checked whether diss_max_step_damping clamps the "
              "dissipation rates");
  }

  if (buffer->ninvalid > 0) {
    CCTK_VERROR("KO dissipation: %.0f points on this process have a "
                "degenerate coordinate map (non-finite or non-positive local "
                "physical spacing h_min = min_d h_d / |grad xi^d|)",
                double(buffer->ninvalid));
  }

  if (buffer->cart_max_step_damping > 0) {
    CCTK_VWARN(CCTK_WARN_ALERT,
               "KO dissipation is clamped on a Cartesian patch on this "
               "process: dt * cart_diss_eps * sum_d 1/h_d = %g exceeds "
               "diss_max_step_damping = %g, so cart_diss_eps is reduced there "
               "by a factor %g. Lower cart_diss_eps or dt to remove the clamp.",
               double(buffer->cart_max_step_damping),
               double(diss_max_step_damping),
               double(buffer->cart_max_step_damping / diss_max_step_damping));
  }

  if (buffer->nclamped > 0) {
    CCTK_VWARN(CCTK_WARN_ALERT,
               "KO dissipation is clamped at %.0f curvilinear-patch points on "
               "this process: dt * sum_d R_d reaches %g, above "
               "diss_max_step_damping = %g, so the rates there are reduced by "
               "up to a factor %g. Smallest local physical spacing among them: "
               "h_min = %g. Lower curv_diss_eps_a,b,c or dt to remove the "
               "clamp.",
               double(buffer->nclamped), double(buffer->max_step_damping),
               double(diss_max_step_damping),
               double(buffer->max_step_damping / diss_max_step_damping),
               double(buffer->min_h_min));
  }
}

} // namespace CapyrX::WaveToy
