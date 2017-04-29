#ifndef CGL_COLLADA_GRIDINFO_H
#define CGL_COLLADA_GRIDINFO_H

#include "collada_info.h"
#include "CGL/spectrum.h"

namespace CGL { namespace Collada {

struct GridInfo : Instance {
  Spectrum sigma_a;
  Spectrum sigma_s;
  float max_density;
  int x;						 ///< x
  int y;						 ///< y
  int z;						 ///< z
  std::vector<float> densities;  ///< density array (may not need)

}; // struct Grid

std::ostream& operator<<(std::ostream& os, const GridInfo& sphere);

} // namespace Collada
} // namespace CGL

#endif // CGL_COLLADA_GRIDINFO_H
