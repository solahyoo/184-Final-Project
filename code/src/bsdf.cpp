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
  return reflectance / PI;
}

Spectrum DiffuseBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // Part 3, Task 1:
  // This function takes in only wo and provides pointers for wi and pdf,
  // which should be assigned by this function.
  // After sampling a value for wi, it returns the evaluation of the BSDF
  // at (wo, *wi).
  *wi = sampler.get_sample(pdf);
  return reflectance / PI;
}


// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO: 3-2 Part 1 Task 2
  // Implement MirrorBSDF
  *pdf = 1.0;
  reflect(wo, wi);
  return reflectance / abs_cos_theta(*wi);
}

// Microfacet BSDF //

Spectrum MicrofacetBSDF::F(const Vector3D& wi) {
    // TODO: proj3-2, part 2
    // Compute Fresnel term for reflection on dielectric-conductor interface.
    // You will need both eta and etaK, both of which are Spectrum.

    double costheta = cos(getTheta(wi));

    Spectrum Rs = (eta * eta + k * k - eta * costheta * 2.0 + Spectrum(costheta * costheta)) /
                  (eta * eta + k * k + eta * costheta * 2.0 + Spectrum(costheta * costheta));

    Spectrum Rp = ((eta * eta + k * k) * costheta * costheta - eta * costheta * 2.0 + Spectrum(1.0)) /
                  ((eta * eta + k * k) * costheta * costheta + eta * costheta * 2.0 + Spectrum(1.0));
    return (Rs + Rp) / 2.0;
}

double MicrofacetBSDF::G(const Vector3D& wo, const Vector3D& wi) {
    return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D& h) {
    // TODO: proj3-2, part 2
    // Compute Beckmann normal distribution function (NDF) here.
    // You will need the roughness alpha.
    double theta = getTheta(h);
    double power = -pow(tan(theta), 2) / pow(alpha, 2);
    double denom = PI * pow(alpha, 2) * pow(cos(theta), 4);
    return exp(power) / denom;

//    return exp(-std::pow(tan(h.z), 2) / std::pow(alpha, 2)) / (PI * std::pow(alpha, 2) * pow(cos(h.z), 4));
//    return std::pow(cos_theta(h), 100.0);
}

Spectrum MicrofacetBSDF::f(const Vector3D& wo, const Vector3D& wi) {
    // TODO: proj3-2, part 2
    // Implement microfacet model here.
    if (wo.z <= 0 || wi.z <= 0) return Spectrum();

    Vector3D h = wo + wi;
    h.normalize();

   return (F(wi) * G(wo, wi) * D(h)) /
           (4.0 * wo.z * wi.z);
}

Spectrum MicrofacetBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
    // TODO: proj3-2, part 2
    // *Importance* sample Beckmann normal distribution function (NDF) here.
    // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
    //       and return the sampled BRDF value.
    Vector2D sample = sampler.get_sample();
    double theta_h = atan(sqrt(-pow(alpha, 2) * log(1 - sample.x)));
    double phi_h = 2.0 * PI * sample.y;

    Vector3D h = Vector3D(sin(theta_h) * cos(phi_h), sin(theta_h) * sin(phi_h), cos(theta_h));
    h.normalize();

    double pdf_theta = (2.0 * sin(theta_h)) / (pow(alpha, 2) * pow(cos(theta_h), 3)) * exp(-pow(tan(theta_h), 2) / pow(alpha, 2));
    double pdf_phi = 1.0 / (2.0 * PI);

    *wi = Vector3D(-wo + 2 * dot(wo, h) * h);
    if (wi->z < 0) {
        *pdf = 0;
        return Spectrum();
    }

    double pdf_h = pdf_phi * pdf_theta / sin(theta_h);
    *pdf = pdf_h / (4.0 * dot(*wi, h));

    if (*pdf <= 0) {
        *pdf = 0;
        return Spectrum();
    }
     return MicrofacetBSDF::f(wo, *wi);

//    *wi = cosineHemisphereSampler.get_sample(pdf);
//    return MicrofacetBSDF::f(wo, *wi);
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
  double n1, n2;
  if (wo.z >= 0) {
    n1 = 1.0;
    n2 = ior;
  } else {
    n1 = ior;
    n2 = 1.0;
  }

  double eta = n1 / n2;

  bool refraction = refract(wo, wi, ior);

  if (!refraction) { // total internal reflection
    reflect(wo, wi);
    *pdf = 1;
    return reflectance / abs_cos_theta(*wi);
  } else {
    double r_0 = pow(((n1 - n2) / (n1 + n2)), 2);
    // double R = clamp(r_0 + (1.0 - r_0) * pow((1 - abs_cos_theta(wo)), 5), 0, 1);
    double R = r_0 + (1.0 - r_0) * pow((1 - abs_cos_theta(wo)), 5);
    if (coin_flip(R)) {
      reflect(wo, wi);
      *pdf = R;
      return *pdf * reflectance / abs_cos_theta(*wi);
    } else {
      refract(wo, wi, ior);
      *pdf = 1.0 - R;
      return *pdf * transmittance / abs_cos_theta(*wi) / pow(eta, 2);
    }
  }

}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {

  // TODO: 3-2 Part 1 Task 1
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
  *wi = Vector3D(-wo[0],-wo[1],wo[2]);

}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {

  // TODO: 3-2 Part 1 Task 3
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When wo.z is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
  double eta;

  if (wo.z > 0) {
    eta = 1.0 / ior;
  } else {
    eta = ior;
//    std::cout << wo.z << eta << ior << std::endl;
  }

  double z = 1 - pow(eta, 2) * (1 - pow(cos_theta(wo), 2));

  if (z < 0) return false;

  *wi = Vector3D(-eta * wo.x, -eta * wo.y, sqrt(z));
  if (wo.z >= 0) {
    wi->z *= -1.0;
  }

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
