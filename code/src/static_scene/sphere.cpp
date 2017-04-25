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
  double a = dot(r.d, r.d);
  double b = dot(2 * (r.o - o), r.d);
  double c = dot(r.o - o, r.o - o) - r2;

  if (std::pow(b, 2) - 4 * a * c < 0) return false;

  double t_1 = (-b - sqrt(pow(b, 2) - 4 * a * c)) / (2.0 * a);
  double t_2 = (-b + sqrt(pow(b, 2) - 4 * a * c)) / (2.0 * a);

  if ((t_1 < r.min_t && t_2 < r.min_t) || (t_1 > r.max_t && t_2 > r.max_t)) return false;

  t1 = t_1;
  t2 = t_2;

  return true;
}

bool Sphere::intersect(const Ray& r) const {

  // Part 1, Task 4:
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.
//  double *t1;
//  double *t2;
//  bool intersects = test(r, (double &) t1, (double &) t2);

  double t1, t2;
  bool intersects = test(r, t1, t2);

  if (!intersects) return false;
//  if (t1 < r.min_t || t1 > r.max_t) return false;

  if (t1 > r.max_t || t2 < r.min_t) return false;

  if (t1 < r.min_t) {
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
  bool intersects = test(r, t1, t2);

  if (!intersects) return false;
//  if (t1 < r.min_t || t1 > r.max_t) return false;

  if (t1 > r.max_t || t2 < r.min_t) return false;

  if (t1 < r.min_t) {
    r.max_t = t2;
    i->t = t2;
  } else {
    r.max_t = t1;
    i->t = t1;
  }

//  r.max_t = t1;
  i->primitive = this;
  i->bsdf = get_bsdf();

  Vector3D p = r.o + i->t * r.d;
  i->n = normal(p);

  return true;

}

float Sphere::medium_dist(const Ray& r) const  {
  double t1, t2;
  bool intersects = test(r, t1, t2);
  if (!intersects) return -1;

  Vector3D d1 = r.o + t1 * r.d;
  Vector3D d2 = r.o + t2 * r.d;

  return (d2 - d1).norm();
  // return 0.0f;
}

// bool Sphere::intersect(const Ray& r, Intersection* i, float& distance) const {
//   return false;
// }

void Sphere::draw(const Color& c) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color& c) const {
    //Misc::draw_sphere_opengl(o, r, c);
}


} // namespace StaticScene
} // namespace CGL
