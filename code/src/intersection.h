#ifndef CGL_INTERSECT_H
#define CGL_INTERSECT_H

#include <vector>

#include "CGL/vector3D.h"
#include "CGL/spectrum.h"
#include "CGL/misc.h"

#include "bsdf.h"

namespace CGL { namespace StaticScene {

class Primitive;
class Grid;

/**
 * A record of an intersection point which includes the time of intersection
 * and other information needed for shading
 */
struct Intersection {

  Intersection() : t (INF_D), primitive(NULL), bsdf(NULL), grid(NULL), is_medium(false) { }

  // constructor for intersections with participating media (with grid)
  Intersection(double t, Grid* grid, const Vector3D &wo) : t(t), grid(grid), wo(wo) {
    primitive = NULL;
    bsdf = NULL;
    n = Vector3D(0, 0, 1);
    is_medium = true;
  }

  double t;    ///< time of intersection

  const Primitive* primitive;  ///< the primitive intersected

  Vector3D n;  ///< normal at point of intersection

  BSDF* bsdf; ///< BSDF of the surface at point of intersection

  Grid* grid; ///< grid of participating media - not null if intersection is with medium

  bool is_medium; ///< true if intersection is with medium

  Vector3D wo;

  // More to follow.
};

} // namespace StaticScene
} // namespace CGL

#endif // CGL_INTERSECT_H
