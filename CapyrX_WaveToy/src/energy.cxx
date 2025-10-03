#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>

namespace CapyrX::WaveToy {

extern "C" void CapyrX_WaveToy_Energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_WaveToy_Energy;
  DECLARE_CCTK_PARAMETERS;

  using namespace Loop;

  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones, [=] CCTK_HOST CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE {
                              const auto l_Pi{Pi(p.I)};
                              const auto l_Dx{Dx(p.I)};
                              const auto l_Dy{Dy(p.I)};
                              const auto l_Dz{Dz(p.I)};

                              energy(p.I) = 0.5 * (l_Pi * l_Pi + l_Dx * l_Dx +
                                                   l_Dy * l_Dy + l_Dz * l_Dz);
                            });
}

} // namespace CapyrX::WaveToy
