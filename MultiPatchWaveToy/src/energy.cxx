// clang-format off
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
// clang-format on

#include "standing_wave.hxx"

namespace MultiPatchWaveToy {

using namespace Arith;

extern "C" void MultiPatchWaveToy_Energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_MultiPatchWaveToy_Energy;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int<0, 0, 0>(grid.nghostzones,
                         [=] CCTK_DEVICE(const Loop::PointDesc &p)
                             CCTK_ATTRIBUTE_ALWAYS_INLINE {
                               // TODO: This contains two nans. Investigate.
                               // const auto l_Pi{Pi(p.I)};
                               // const auto l_Dx{Dx(p.I)};
                               // const auto l_Dy{Dy(p.I)};
                               // const auto l_Dz{Dz(p.I)};

                               // energy(p.I) = 0.5 * (l_Pi * l_Pi + l_Dx * l_Dx
                               // + l_Dy * l_Dy + l_Dz * l_Dz);

                               energy(p.I) = 0.0;
                             });
}

} // namespace MultiPatchWaveToy
