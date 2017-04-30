#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CGL {

bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

  // Part 2, Task 2:
  // Implement ray - bounding box intersection test
  // If the ray intersected the bouding box within the range given by
  // t0, t1, update t0 and t1 with the new intersection times.

  double t_x1 = (min.x - r.o.x) / r.d.x;
  double t_x2 = (max.x - r.o.x) / r.d.x;
  double t_y1 = (min.y - r.o.y) / r.d.y;
  double t_y2 = (max.y - r.o.y) / r.d.y;
  double t_z1 = (min.z - r.o.z) / r.d.z;
  double t_z2 = (max.z - r.o.z) / r.d.z;

  double t_xmin, t_xmax, t_ymin, t_ymax, t_zmin, t_zmax;
  if (t_x1 < t_x2) {
    t_xmin = t_x1;
    t_xmax = t_x2;
  } else {
    t_xmin = t_x2;
    t_xmax = t_x1;
  }

  if (t_y1 < t_y2) {
    t_ymin = t_y1;
    t_ymax = t_y2;
  } else {
    t_ymin = t_y2;
    t_ymax = t_y1;
  }

  if (t_z1 < t_z2) {
    t_zmin = t_z1;
    t_zmax = t_z2;
  } else {
    t_zmin = t_z2;
    t_zmax = t_z1;
  }

  double t_min = std::max(t_xmin, std::max(t_ymin, t_zmin));
  double t_max = std::min(t_xmax, std::min(t_ymax, t_zmax));

  if ((t_min < r.min_t && t_max < r.min_t) ||  (t_min > r.max_t && t_max > r.max_t)) return false;
  if (t_min <= t_max) {
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
