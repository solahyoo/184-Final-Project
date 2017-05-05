#include "grid.h"
#include "CGL/vector3D.h"

#include <iostream>

using namespace std;
namespace CGL { namespace StaticScene {

  float generate_rand() {
    return float(rand()) / RAND_MAX;
  }

  double lerp(double x, double v0, double v1) {
    return (1 - x) * v0 + x * v1;
  }

  double Grid::trilerp_density(const Vector3D& v) const {
    // Compute coordinates and offsets for v
    Vector3D samples = Vector3D(v.x * x - 0.5, v.y * y - 0.5, v.z * z - 0.5);
    Vector3D vi = Vector3D(floor(samples.x), floor(samples.y), floor(samples.z));
    Vector3D d = samples - vi;

    // Trilinearly interpolate density values to compute local density
    double d00 = lerp(d.x, D(vi), D(vi + Vector3D(1, 0, 0)));
    double d10 = lerp(d.x, D(vi + Vector3D(0, 1, 0)), D(vi + Vector3D(1, 1, 0)));
    double d01 = lerp(d.x, D(vi + Vector3D(0, 0, 1)), D(vi + Vector3D(1, 0, 1)));
    double d11 = lerp(d.x, D(vi + Vector3D(0, 1, 1)), D(vi + Vector3D(1, 1, 1)));
    double d0 = lerp(d.y, d00, d10);
    double d1 = lerp(d.y, d01, d11);
    return lerp(d.z, d0, d1);
  }

  Spectrum Grid::sample(const Ray& r, Intersection *i) {
    Ray mray = Ray(r.o, r.d.unit());
    if (mray.max_t != INF_D)
      mray.max_t = r.max_t * r.d.norm();
    mray.o = w2g * mray.o;
    mray.d = w2g * mray.d;
    double tmin, tmax;
    BBox b = get_bbox();
    if (!b.intersect(mray, tmin, tmax))
      return Spectrum(1, 1, 1);
    double t = tmin;
    while (true) {
      float random = generate_rand();
      t -= std::log(1 - random) / (max_density * sigma_t);
      if (t >= tmax)
        break;
      if (trilerp_density(mray.o + mray.d * t) / max_density > generate_rand()) {
        *i = Intersection(tmin, this);
        return sigma_s / sigma_t;
      }
    }
    return Spectrum(1, 1, 1);
  }

  Spectrum Grid::transmittance(const Ray& r) const {
    Ray mray = Ray(r.o, r.d.unit());
    if (mray.max_t != INF_D)
      mray.max_t = r.max_t * r.d.norm();
    mray.o = w2g * mray.o;
    mray.d = w2g * mray.d;
    double tmin, tmax;
    BBox b = get_bbox();
    if (!b.intersect(mray, tmin, tmax))
      return Spectrum(1, 1, 1);
    double tr = 1;
    double t = tmin;
    while (true) {
      float random = generate_rand();
      t -= std::log(1 - random) / (max_density * sigma_t);
      if (t >= tmax)
        break;
      double d = trilerp_density(mray.o + mray.d * t);
      tr *= 1 - std::max(0.0, d / max_density);

      // when transmittance gets low, start applying Russian roulette to terminate sampling
      if (tr < .1) {
        double a = std::max(0.05, 1 - tr);
        if (generate_rand() < a)
          return Spectrum(0.05, 0.05, 0.05);
        tr /= 1 - a;
      }
    }
    return Spectrum(tr, tr, tr);
  }

  double Grid::p(const Vector3D& wo, const Vector3D& wi) {
    return phaseHG(dot(wo, wi));
  }
  double Grid::sample_p(const Vector3D &wo, Vector3D *wi, const Vector2D &u) {
    double cosTheta;
    if (std::abs(g) < 1e-3)
      cosTheta = 1 - 2 * u.x;
    else {
      double sqrTerm = (1 - g * g) / (1 - g + 2 * g * u.x);
      cosTheta = (1 + g * g - sqrTerm * sqrTerm) / (2 * g);
    }

    // Compute direction _wi_ for Henyey--Greenstein sample
    double sinTheta = std::sqrt(std::max(0.0, 1 - cosTheta * cosTheta));
    double phi = 2 * PI * u.y;
    Vector3D v1, v2;
    if (std::abs(wo.x) > std::abs(wo.y))
      v1 = Vector3D(-wo.z, 0, wo.x) / std::sqrt(wo.x * wo.x + wo.z * wo.z);
    else
      v1 = Vector3D(0, wo.z, -wo.y) / std::sqrt(wo.y * wo.y + wo.z * wo.z);
    v2 = cross(wo, v1);
    // 1 / 4pi
    // spherical direction
    *wi = sinTheta * cos(phi) * v1 + sinTheta * sin(phi) * v2 + cosTheta * -wo;
    return phaseHG(-cosTheta);
  }


}}
