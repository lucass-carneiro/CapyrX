#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>
#include <global_derivatives.hxx>

#include "local_derivatives.hxx"

namespace MultiPatchWaveToy2 {

// u(t,x,y,z) =
//   A cos(2 pi omega t) sin(2 pi kx x) sin(2 pi ky y) sin(2 pi kz z)
template <typename T>
constexpr void standing_wave(const T A, const T kx, const T ky, const T kz,
                             const T t, const T x, const T y, const T z, T &u,
                             T &rho) {
  using std::acos, std::cos, std::pow, std::sin, std::sqrt;

  const T pi = acos(-T(1));
  const T omega = sqrt(pow(kx, 2) + pow(ky, 2) + pow(kz, 2));

  u = A * cos(2 * pi * omega * t) * cos(2 * pi * kx * x) *
      cos(2 * pi * ky * y) * cos(2 * pi * kz * z);
  rho = A * (-2 * pi * omega) * sin(2 * pi * omega * t) * cos(2 * pi * kx * x) *
        cos(2 * pi * ky * y) * cos(2 * pi * kz * z);
}

// u(t,r) = (f(t-r) - f(t+r)) / r
// f(v) = A exp(-1/2 (r/W)^2)
template <typename T>
constexpr void gaussian(const T A, const T W, const T t, const T x, const T y,
                        const T z, T &u, T &rho) {
  using std::exp, std::pow, std::sqrt;

  const T r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

  const auto f = [&](const T v) {
    return A * exp(-pow(v, 2) / (2 * pow(W, 2)));
  };

  if (r < sqrt(std::numeric_limits<T>::epsilon())) {
    // L'Hôpital
    u = 2 / pow(W, 2) * f(t) * t;
    rho = -2 / pow(W, 4) * f(t) * (pow(t, 2) - pow(W, 2));
  } else {
    u = (f(t - r) - f(t + r)) / r;
    rho = -(f(t - r) * (t - r) - f(t + r) * (t + r)) / (pow(W, 2) * r);
  }
}

extern "C" void MultiPatchWaveToy2_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_MultiPatchWaveToy2_Initial;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto t{cctk_time};
          const auto x{vcoordx(p.I)};
          const auto y{vcoordy(p.I)};
          const auto z{vcoordz(p.I)};

          standing_wave(amplitude, wave_kx, wave_ky, wave_kz, t, x, y, z,
                        u(p.I), rho(p.I));
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto t{cctk_time};
          const auto x{vcoordx(p.I)};
          const auto y{vcoordy(p.I)};
          const auto z{vcoordz(p.I)};

          gaussian(amplitude, gaussian_width, t, x, y, z, u(p.I), rho(p.I));
        });
  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

extern "C" void MultiPatchWaveToy2_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_MultiPatchWaveToy2_RHS;
  DECLARE_CCTK_PARAMETERS;

  using namespace MultiPatch::GlobalDerivatives;

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const LocalFirstDerivatives ldu{c4o_1_0_0(p, u), c4o_0_1_0(p, u),
                                        c4o_0_0_1(p, u)};

        const LocalSecondDerivatives ld2u{c4o_2_0_0(p, u), c4o_1_1_0(p, u),
                                          c4o_1_0_1(p, u), c4o_0_2_0(p, u),
                                          c4o_0_1_1(p, u), c4o_0_0_2(p, u)};

        const Jacobians jac{VERTEX_JACOBIANS_FIRST(p)};
        const JacobianDerivatives djac{VERTEX_JACOBIANS_SECOND(p)};

        const auto d2u{project_second(ldu, ld2u, jac, djac)};

        u_rhs(p.I) = rho(p.I);
        rho_rhs(p.I) = d2u.dxdx + d2u.dydy + d2u.dzdz;
      });
}

extern "C" void MultiPatchWaveToy2_Sync(CCTK_ARGUMENTS) {
  // Do nothing
}

extern "C" void MultiPatchWaveToy2_Energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_MultiPatchWaveToy2_Energy;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const Loop::PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      // TODO
                                      energy(p.I) = 0.0;
                                    });
}

extern "C" void MultiPatchWaveToy2_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_MultiPatchWaveToy2_Error;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          using std::fabs;

          const auto t{cctk_time};
          const auto x{vcoordx(p.I)};
          const auto y{vcoordy(p.I)};
          const auto z{vcoordz(p.I)};

          CCTK_REAL u0, rho0;
          standing_wave(amplitude, wave_kx, wave_ky, wave_kz, t, x, y, z, u0,
                        rho0);

          u_error(p.I) = fabs(u(p.I) - u0);
          rho_error(p.I) = fabs(rho(p.I) - rho0);
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          using std::fabs;

          const auto t{cctk_time};
          const auto x{vcoordx(p.I)};
          const auto y{vcoordy(p.I)};
          const auto z{vcoordz(p.I)};

          CCTK_REAL u0, rho0;
          gaussian(amplitude, gaussian_width, t, x, y, z, u0, rho0);

          u_error(p.I) = fabs(u(p.I) - u0);
          rho_error(p.I) = fabs(rho(p.I) - rho0);
        });
  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

} // namespace MultiPatchWaveToy2
