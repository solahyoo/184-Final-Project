#include "grid.h"


using namespace std;
namespace CGL { namespace StaticScene {

  Spectrum Grid::transmittance(const Ray& r) const {
    return Spectrum();
  }

  Spectrum Grid::sample(const Ray& r) {
    Ray med_ray = Ray(r.o, r.d.unit())
    med_ray.max_t = r.max_t * r.d.norm();
    float *tmin, *tmax;
    BBox b = get_bbox();
    if (!b.Intersect(med_ray, &tmin, &tmax))
      return Spectrum();
    float t = tmin;
    while (true) {
      // translating pbrt code to ours
      float random = float(rand()) / RAND_MAX;
      t -= std::log(1 - random) / (max_density * sigma_t);
      if (t >= tmax)
        break;
      // if () {
      //   // create phase function
      //   return sigma_s / sigma_t;
      // }
    }
    return Spectrum();
  }

}}
