#include "grid.h"
#include "CGL/vector3D.h"


using namespace std;
namespace CGL { namespace StaticScene {

  float generate_rand() {
    return float(rand()) / RAND_MAX;
  }

  float lerp(float x, float v0, float v1) {
    return (1 - x) * v0 + x * v1;
  }

  float Grid::trilerp_density(const Vector3D& v) const {
    // Compute coordinates and offsets for v
    Vector3D samples = Vector3D(v.x * x - 0.5, v.y * y - 0.5, v.z * z - 0.5);
    Vector3D vi = Vector3D(floor(samples.x), floor(samples.y), floor(samples.z));
    Vector3D d = samples - vi;

    // Trilinearly interpolate density values to compute local density
    float d00 = lerp(d.x, D(vi), D(vi + Vector3D(1, 0, 0)));
    float d10 = lerp(d.x, D(vi + Vector3D(0, 1, 0)), D(vi + Vector3D(1, 1, 0)));
    float d01 = lerp(d.x, D(vi + Vector3D(0, 0, 1)), D(vi + Vector3D(1, 0, 1)));
    float d11 = lerp(d.x, D(vi + Vector3D(0, 1, 1)), D(vi + Vector3D(1, 1, 1)));
    float d0 = lerp(d.y, d00, d10);
    float d1 = lerp(d.y, d01, d11);
    return lerp(d.z, d0, d1);
  }

  Spectrum Grid::sample(const Ray& r) {
    Ray mray = Ray(r.o, r.d.unit());
    mray.max_t = r.max_t * r.d.norm();
    double tmin, tmax;
    BBox b = get_bbox();
    if (!b.intersect(mray, tmin, tmax))
      return Spectrum(1, 1, 1);
    float t = tmin;
    while (true) {
      // translating pbrt code to ours
      float random = generate_rand();
      t -= std::log(1 - random) / (max_density * sigma_t);
      if (t >= tmax)
        break;
      if (trilerp_density(mray.o + mray.d * t) / max_density > generate_rand()) {
        // create phase function
        return sigma_s / sigma_t;
      }
    }
    return Spectrum(1, 1, 1);
  }

  Spectrum Grid::transmittance(const Ray& r) const {
    Ray mray = Ray(r.o, r.d.unit());
    mray.max_t = r.max_t * r.d.norm();
    double tmin, tmax;
    BBox b = get_bbox();
    if (!b.intersect(mray, tmin, tmax))
      return Spectrum(1, 1, 1);
    float tr = 1;
    float t = tmin;
    while (true) {
      float random = generate_rand();
      t -= std::log(1 - random) / (max_density * sigma_t);
      if (t >= tmax)
        break;
      float density = trilerp_density(mray.o + mray.d * t);
      tr *= 1 - std::max(0.0f, density / max_density);

      // when trnasmittance gets low, start applying Russian roulette to terminate sampling
      if (tr < .1) {
        float a = std::max(0.05f, 1 - tr);
        if (generate_rand() < a)
          return Spectrum();
        tr /= 1 - a;
      }
    }
    return Spectrum(tr, tr, tr);
  }

}}
