#include "sphere.h"

#include <cmath>

#include  "../bsdf.h"
#include "../misc/sphere_drawing.h"

namespace CGL { namespace StaticScene {

bool Sphere::test(const Ray& r, double& t1, double& t2) const {

  // Part 1, Task 4:
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.
  Vector3D roo = r.o - o;
  double t;
  double a = dot(r.d, r.d);
  double b = dot(2 * roo, r.d);
  double c = dot(roo, roo) - r2;
  double determinant = b * b - 4 * a * c;
  if (determinant < 0) {
    return false;
  }
  t1 = (-b - sqrt(determinant)) / (2 * a);
  t2 = (-b + sqrt(determinant)) / (2 * a);
  if (t1 > t2) {
    double temp = t1;
    t1 = t2;
    t2 = temp;
  }
  return true;
}

bool Sphere::intersect(const Ray& r) const {

  // Part 1, Task 4:
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.
  double t1, t2;
  bool intersection = test(r, t1, t2);
  if (!intersection || t1 > t2)
    return false;
  if (t1 > r.max_t || t2 < r.min_t)
    return false;

  if (t1 < r.min_t && t2 <= r.max_t) {
    r.max_t = t2;
  } else {
    r.max_t = t1;
  }
  return true;
}

bool Sphere::intersect(const Ray& r, Intersection *i) const {

  // Part 1, Task 4:
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.
  double t1, t2;
  bool intersection = test(r, t1, t2);
  // check that t1 and t2 are valid
  if (!intersection || t1 > t2)
    return false;
  // if (t1 > t2 || !intersection)
  // can now assume t1 <= t2
  if (t1 > r.max_t || t2 < r.min_t)
    return false;

  if (t1 < r.min_t && t2 <= r.max_t) {
    i->t = t2;
    r.max_t = t2;
  } else {
    i->t = t1;
    r.max_t = t1;
  }
  i->n = (r.o + i->t * r.d - o).unit();
  i->primitive = this;
  i->bsdf = get_bsdf();
  return true;
}

void Sphere::draw(const Color& c) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color& c) const {
    //Misc::draw_sphere_opengl(o, r, c);
}


} // namespace StaticScene
} // namespace CGL
