#ifndef CGL_GRID_H
#define CGL_GRID_H

#include "CGL/CGL.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"
#include "CGL/matrix4x4.h"
#include "CGL/spectrum.h"
#include "collada/grid_info.h"
#include "math.h"
#include "ray.h"
#include "bbox.h"
#include "intersection.h"
#include <vector>
#include <cstdio>

using namespace std;

namespace CGL { namespace StaticScene {

/**
 * Class for grid for participating media
 * Want to represent grid of densities.
 */
class Grid {
 public:
   /**
    * Constructor.
    */
  Grid(Collada::GridInfo& grid_info, const Matrix4x4& transform) {
    sigma_a = grid_info.sigma_a;
    sigma_s = grid_info.sigma_s;
    max_density = grid_info.max_density;
    x = grid_info.x;
    y = grid_info.y;
    z = grid_info.z;
    density = grid_info.densities;
    sigma_t = (sigma_a + sigma_s).r;

    g = .5; // change later
  }
  Grid(const Spectrum &sigma_a, const Spectrum &sigma_s, float max_density, int x, int y, int z, vector<float> d)
    : sigma_a(sigma_a),
      sigma_s(sigma_s),
      max_density(max_density),
      x(x), y(y), z(z),
      density(d) {

      // memcpy((float *)density.get(), d, sizeof(float) * x * y * z);
      sigma_t = (sigma_a + sigma_s).r;

      g = .5; // change later
    }

  BBox get_bbox() const {
    return BBox(Vector3D(), Vector3D(1, 1, 1));
  }

  float D(const Vector3D &v) const {
    bool inside_exclusive = v.x >= 0 && v.x < x && v.y >= 0 && v.y < y &&
                            v.z >= 0 && v.z < z;
    if (!inside_exclusive)
      return 0;
    return density[(v.z * y + v.y) * x + v.x];
  }

  float intersection(const Ray& r) const;

  float trilerp_density(const Vector3D& v) const;

  Spectrum sample(const Ray& r, Intersection *i);

  Spectrum transmittance(const Ray& r) const;

  // For phase function - Henyey Greenstein phase function
  float phaseHG(float cosTheta) {
    float denom = 1 + g * g + 2 * g * cosTheta;
    return (1 - g * g) / (denom * std::sqrt(denom) * 4 * PI);
  }
  float p(const Vector3D& wo, const Vector3D& wi);
  float sample_p(const Vector3D &wo, Vector3D *wi, const Vector2D &u);



private:
 // dimensions of the grid
 int x, y, z;
 float max_density;

 // coefficient for Henyey Greenstein phase function
 float g;

 // array of densities for the grid
 vector<float> density;

 // absorption coefficient
 Spectrum sigma_a;
 // scattering coefficient
 Spectrum sigma_s;
 float sigma_t;

 // world to grid rotation matrix
 Matrix3x3 w2g;

 };

}} // namespace CGL

#endif // CGL_GRID_H