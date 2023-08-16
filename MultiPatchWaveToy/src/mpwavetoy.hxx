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

constexpr const int dim = 3;

extern "C" {

void MultiPatchWaveToy_Initial(CCTK_ARGUMENTS);
void MultiPatchWaveToy_RHS(CCTK_ARGUMENTS);
void MultiPatchWaveToy_Energy(CCTK_ARGUMENTS);
void MultiPatchWaveToy_Error(CCTK_ARGUMENTS);
}

} // namespace MultiPatchWaveToy

#endif