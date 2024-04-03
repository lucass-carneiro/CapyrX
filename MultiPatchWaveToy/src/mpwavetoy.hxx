#ifndef CAPYRX_MULTIPATCH_WAVE_TOY_HPP
#define CAPYRX_MULTIPATCH_WAVE_TOY_HPP

// clang-format off
#include <loop_device.hxx>

#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
//clang-format on

namespace MultiPatchWaveToy {

constexpr int dim{3};

/*
 * u(t,x,y,z) = A cos(2 pi omega t) sin(2 pi kx x) sin(2 pi ky y) sin(2 pi kz z)
 */
void standing_wave(const CCTK_REAL A, const CCTK_REAL kx, const CCTK_REAL ky, const CCTK_REAL kz,
                             const CCTK_REAL t, const CCTK_REAL x, const CCTK_REAL y, const CCTK_REAL z, CCTK_REAL &u,
                             CCTK_REAL &rho) noexcept;
                             
/*
 * u(t,r) = (f(t-r) - f(t+r)) / r
 * f(v) = A exp(-1/2 (r/W)^2)
 */
void gaussian(const CCTK_REAL A, const CCTK_REAL W, const CCTK_REAL t, const CCTK_REAL x, const CCTK_REAL y,
                        const CCTK_REAL z, CCTK_REAL &u, CCTK_REAL &rho) noexcept;

extern "C" {

void MultiPatchWaveToy_Initial(CCTK_ARGUMENTS);
void MultiPatchWaveToy_RHS(CCTK_ARGUMENTS);
void MultiPatchWaveToy_Error(CCTK_ARGUMENTS);

}

} // namespace MultiPatchWaveToy

#endif