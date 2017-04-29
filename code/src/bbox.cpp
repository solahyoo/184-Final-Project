#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CGL {

/* Divides each element of vect1 by corresponding element of vect2 */
Vector3D div_each(const Vector3D v1, const Vector3D v2) {
  return Vector3D(v1[0] / v2[0], v1[1] / v2[1], v1[2] / v2[2]);
}

/* Returns the min of v1[index] or v2[index] if find_max = 0
   Returns the max if find_max = 1
*/
double find_min_max(Vector3D v1, Vector3D v2, int index, int find_max) {
  if (find_max)
    return std::max(v1[index], v2[index]);
  return std::min(v1[index], v2[index]);
}

/* Returns the max of elements in v1 if find_max == 1
   Else, returns the min of elements
*/
double vect_min_max(Vector3D v1, int find_max) {
  if (find_max)
    return std::max(v1.x, std::max(v1.y, v1.z));
  return std::min(v1.x, std::min(v1.y, v1.z));
}

bool BBox::intersect(const Ray& r, double& t0, double& t1) const {
  // Part 2, Task 2:
  // Implement ray - bounding box intersection test
  // If the ray intersected the bouding box within the range given by
  // t0, t1, update t0 and t1 with the new intersection times.
  Vector3D v1 = div_each((min - r.o), r.d);
  Vector3D v2 = div_each((max - r.o), r.d);
  Vector3D vmin = Vector3D(find_min_max(v1, v2, 0, 0), find_min_max(v1, v2, 1, 0), find_min_max(v1, v2, 2, 0));
  Vector3D vmax = Vector3D(find_min_max(v1, v2, 0, 1), find_min_max(v1, v2, 1, 1), find_min_max(v1, v2, 2, 1));
  double tmin = vect_min_max(vmin, 1); // get max of vmin
  double tmax = vect_min_max(vmax, 0); // get min of vmax

  // check if t0 or t1 are valid intersection points
  if ((tmin <= r.max_t || tmax >= r.min_t) && tmin <= tmax) {
    t0 = tmin;
    t1 = tmax;
    return true;
  }
  return false;
}

void BBox::draw(Color c) const {

  glColor4f(c.r, c.g, c.b, c.a);

	// top
	glBegin(GL_LINE_STRIP);
	glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(max.x, max.y, max.z);
	glEnd();

	// bottom
	glBegin(GL_LINE_STRIP);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, min.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, min.y, min.z);
	glEnd();

	// side
	glBegin(GL_LINES);
	glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, min.y, max.z);
	glVertex3d(max.x, max.y, min.z);
  glVertex3d(max.x, min.y, min.z);
	glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, min.y, min.z);
	glVertex3d(min.x, max.y, max.z);
  glVertex3d(min.x, min.y, max.z);
	glEnd();

}

std::ostream& operator<<(std::ostream& os, const BBox& b) {
  return os << "BBOX(" << b.min << ", " << b.max << ")";
}

} // namespace CGL
