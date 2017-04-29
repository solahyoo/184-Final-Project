#include "grid_info.h"

namespace CGL { namespace Collada {

std::ostream& operator<<(std::ostream& os, const GridInfo& grid) {
  return os << "GridInfo: " << grid.name << " (id: " << grid.id << ")"
            << " ["
            << " sigma_a="  << grid.sigma_a
            << " sigma_s="  << grid.sigma_s
            << " x="  << grid.x
            << " y="  << grid.y
            << " z="  << grid.z
            << " num_densities="  << grid.densities.size()
            << " ]";
}

} // namespace Collada
} // namespace CGL
