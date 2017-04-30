#include "triangle.h"

#include "CGL/CGL.h"
#include "GL/glew.h"

namespace CGL { namespace StaticScene {

Triangle::Triangle(const Mesh* mesh, size_t v1, size_t v2, size_t v3) :
    mesh(mesh), v1(v1), v2(v2), v3(v3) { }

BBox Triangle::get_bbox() const {

  Vector3D p1(mesh->positions[v1]), p2(mesh->positions[v2]), p3(mesh->positions[v3]);
  BBox bb(p1);
  bb.expand(p2);
  bb.expand(p3);
  return bb;

}

Vector3D get_intersection(const Ray& r, Vector3D p1, Vector3D p2, Vector3D p3) {
  Vector3D E1 = p2 - p1;
  Vector3D E2 = p3 - p1;
  Vector3D S = r.o - p1;
  Vector3D S1 = cross(r.d, E2);
  Vector3D S2 = cross(S, E1);
  Vector3D result = (1.0 / dot(S1, E1)) * Vector3D(dot(S2, E2), dot(S1, S), dot(S2, r.d));
  return result;
}

bool Triangle::intersect(const Ray& r) const {

  // Part 1, Task 3: implement ray-triangle intersection
  Vector3D p1(mesh->positions[v1]), p2(mesh->positions[v2]), p3(mesh->positions[v3]);
  Vector3D result = get_intersection(r, p1, p2, p3);
  double t = result[0];
  if (t < r.min_t || t > r.max_t)
    return false;
  bool intersect = t >= 0 && result[1] >= 0 && result[1] <= 1 &&
                    result[2] >= 0 && result[2] <= 1 &&
                    (1 - result[1] - result[2]) >= 0 && (1 - result[1] - result[2]) <= 1;
  if (!intersect) return false;
  r.max_t = t;
  return true;
}

bool Triangle::intersect(const Ray& r, Intersection *isect) const {

  // Part 1, Task 3:
  // implement ray-triangle intersection. When an intersection takes
  // place, the Intersection data should be updated accordingly
  Vector3D p1(mesh->positions[v1]), p2(mesh->positions[v2]), p3(mesh->positions[v3]);
  Vector3D n1(mesh->normals[v1]), n2(mesh->normals[v2]), n3(mesh->normals[v3]);
  Vector3D result = get_intersection(r, p1, p2, p3);
  double t = result[0];
  if (t < r.min_t || t > r.max_t)
    return false;
  // get barycentric coordinates
  double gamma = 1 - result[1] - result[2];
  if (result[1] < 0 || result[2] < 0 || gamma < 0)
    return false;
  r.max_t = t;
  Vector3D barycentric = gamma * n1 + result[1] * n2 + result[2] * n3;
  isect->n = barycentric;
  isect->t = t;
  isect->primitive = this;
  isect->bsdf = get_bsdf();
  return true;
}

void Triangle::draw(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_TRIANGLES);
  glVertex3d(mesh->positions[v1].x,
             mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x,
             mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x,
             mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}

void Triangle::drawOutline(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_LINE_LOOP);
  glVertex3d(mesh->positions[v1].x,
             mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x,
             mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x,
             mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}



} // namespace StaticScene
} // namespace CGL
