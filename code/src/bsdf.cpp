#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>

using std::min;
using std::max;
using std::swap;

namespace CGL {

void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {

    Vector3D z = Vector3D(n.x, n.y, n.z);
    Vector3D h = z;
    if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
    else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
    else h.z = 1.0;

    z.normalize();
    Vector3D y = cross(h, z);
    y.normalize();
    Vector3D x = cross(z, y);
    x.normalize();

    o2w[0] = x;
    o2w[1] = y;
    o2w[2] = z;
}


// Diffuse BSDF //

Spectrum DiffuseBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  // Part 3, Task 1:
  // This function takes in both wo and wi and returns the evaluation of
  // the BSDF for those two directions.
  // wi = incoming light direction
  return reflectance / PI;
}

Spectrum DiffuseBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // Part 3, Task 1:
  // This function takes in only wo and provides pointers for wi and pdf,
  // which should be assigned by this function.
  // After sampling a value for wi, it returns the evaluation of the BSDF
  // at (wo, *wi).
  *wi = sampler.get_sample(pdf);
  return f(wo, *wi);
}


// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO: 3-2 Part 1 Task 2
  // Implement MirrorBSDF
  BSDF::reflect(wo, wi);
  *pdf = 1;
  return reflectance / abs_cos_theta(*wi);
}


// Microfacet BSDF //

Spectrum MicrofacetBSDF::F(const Vector3D& wi) {
    // TODO: proj3-2, part 2
    // Compute Fresnel term for reflection on dielectric-conductor interface.
    // You will need both eta and K, both of which are Spectrum.
    Spectrum eta_k = eta * eta + k * k;
    // Spectrum rs = (eta_k - 2 * eta * wi.z + wi.z * wi.z) / (eta_k + 2 * eta * wi.z + wi.z * wi.z);
    // Spectrum rp = (eta_k * wi.z * wi.z - 2 * eta * wi.z + 1) / (eta_k * wi.z * wi.z + 2 * eta * wi.z + 1);
    Spectrum rs = (eta_k - 2 * eta * wi.z + pow(wi.z, 2)) / (eta_k + 2 * eta * wi.z + pow(wi.z, 2));
    Spectrum rp = (eta_k * pow(wi.z, 2) - 2 * eta * wi.z + 1) / (eta_k * pow(wi.z, 2) + 2 * eta * wi.z + 1);
    return (rs + rp) / 2.0;
}

double MicrofacetBSDF::G(const Vector3D& wo, const Vector3D& wi) {
    return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D& h) {
    // TODO: proj3-2, part 2
    // Compute Beckmann normal distribution function (NDF) here.
    // You will need the roughness alpha.
    // theta_h = angle between h and macro surface normal n

    double num = exp(-(1 - pow(h.z, 2)) / (pow(h.z, 2) * pow(alpha, 2)));
    double denom =  (PI * pow(alpha, 2) * pow(h.z, 4));
    return num / denom;
}

Spectrum MicrofacetBSDF::f(const Vector3D& wo, const Vector3D& wi) {
    // TODO: proj3-2, part 2
    // Implement microfacet model here.

    // check if wi and wo are valid. If not, return 0.
    if (wi.z <= 0 || wo.z <= 0) {
      return Spectrum();
    }
    Vector3D h = (wo.unit() + wi.unit()).unit();
    Spectrum s = (F(wi) * G(wo, wi) * D(h)) / (4 * wo.z * wi.z);
    return s;
}

Spectrum MicrofacetBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
    // TODO: proj3-2, part 2
    // *Importance* sample Beckmann normal distribution function (NDF) here.
    // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
    //       and return the sampled BRDF value.

    // *wi = cosineHemisphereSampler.get_sample(pdf);
    // return MicrofacetBSDF::f(wo, *wi);

    Vector2D r = sampler.get_sample();
    double theta_h = atan(sqrt(-pow(alpha, 2) * log(1 - r.x)));
    double phi_h = 2 * PI * r.y;
    Vector3D h = Vector3D(cos(phi_h) * sin(theta_h), sin(phi_h) * sin(theta_h), cos(theta_h));
    *wi = -wo + 2 * dot(wo, h) * h;
    double p_theta_h = 2 * exp(-pow(tan(theta_h) / alpha, 2)) * sin(theta_h) /
                          (pow(alpha, 2) * pow(cos(theta_h), 3));
    double p_phi_h = 1.0 / (2.0 * PI);
    double pdf_h = p_theta_h * p_phi_h / sin(theta_h);
    *pdf = pdf_h / (4 * dot(*wi, h));
    if (wi->z <= 0 || *pdf <= 0) {
      *pdf = 0;
      return Spectrum();
    }
    return MicrofacetBSDF::f(wo, *wi);
}


// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {


  return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO: 3-2 Part 1 Task 4
  // Compute Fresnel coefficient and either reflect or refract based on it.

  // returns false if refraction does not occur due to total internal reflection
  bool refraction = BSDF::refract(wo, wi, ior);
  double eta = ior;
  double n1 = 1.0;
  double n2 = ior;
  if (wo.z > 0) {
    n1 = ior;
    n2 = 1.0;
    eta = 1.0 / ior;
  }
  if (!refraction) {
    // assign reflection of wo to *wi
    BSDF::reflect(wo, wi);
    *pdf = 1;
    return reflectance / abs_cos_theta(*wi);
  } else {
    double r0 = pow((n1 - n2) / (n1 + n2), 2);
    double R = clamp(r0 + (1 - r0) * pow((1 - abs_cos_theta(wo)), 5), 0.0, 1.0);
    if (coin_flip(R)) {
      BSDF::reflect(wo, wi);
      *pdf = R;
      return R * reflectance / abs_cos_theta(*wi);
    } else {
      // refraction of wo to *wi is already assigned in first line
      *pdf = 1.0 - R;
      return (1.0 - R) * transmittance / abs_cos_theta(*wi) / pow(eta, 2);
    }
  }
}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {
  // TODO: 3-2 Part 1 Task 1
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
  *wi = Vector3D(-wo[0], -wo[1], wo[2]);
}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {

  // TODO: 3-2 Part 1 Task 3
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When wo.z is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
  double eta;
  // wo is entering non-air material
  if (wo.z > 0) {
    // ior is new index of refraction of material wi is pointing to
    eta = 1.0 / ior;
  } else { // wo is exiting - when wo starts out inside object with ior > 0
    // ior is old index of refraction of material that wo is inside of
    eta = ior;
  }
  if ((1 - eta * eta * (1 - wo.z * wo.z)) < 0)
    return false;

  double z = sqrt(1 - eta * eta * (1 - wo.z * wo.z));
  if ((wo.z > 0 && z > 0) || (wo.z < 0 && z < 0)) {
    z *= -1.0;
  }
  *wi = Vector3D(-eta * wo.x, -eta * wo.y, z);
  return true;
}

// Emission BSDF //

Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0 / PI;
  *wi  = sampler.get_sample(pdf);
  return Spectrum();
}

} // namespace CGL
