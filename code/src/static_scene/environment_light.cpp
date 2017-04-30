#include "environment_light.h"

#include <algorithm>
#include <iostream>
#include <fstream>

namespace CGL { namespace StaticScene {

EnvironmentLight::EnvironmentLight(const HDRImageBuffer* envMap)
    : envMap(envMap) {
    	init();
}

EnvironmentLight::~EnvironmentLight() {
    delete[] pdf_envmap;
    delete[] conds_y;
    delete[] marginal_y;
}


void EnvironmentLight::init() {
	uint32_t w = envMap->w, h = envMap->h;
    pdf_envmap = new double[w * h];
	conds_y = new double[w * h];
  marginal_y = new double[h];
	double* prob_y = new double[h];

	std::cout << "[PathTracer] Initializing environment light...";

	// TODO 3-2 Part 3 Task 3 Steps 1,2,3
	// Store the environment map pdf to pdf_envmap
	// Store the marginal distribution for y to marginal_y
	// Store the conditional distribution for x given y to conds_y
  double sum = 0;
  for (int j = 0; j < h; ++j) {
		for (int i = 0; i < w; ++i) {
      sum += envMap->data[j * w + i].illum() * sin(PI * j / h);
		}
	}
  for (int j = 0; j < h; ++j) {
    double prob_j = 0;
		for (int i = 0; i < w; ++i) {
      pdf_envmap[j * w + i] = envMap->data[j * w + i].illum() * sin(PI * (j) / h) / sum;
      prob_j += pdf_envmap[j * w + i];
		}
    prob_y[j] = prob_j;
	}
  double marginal_distr = 0;
  double curr_pdf;
  for (int j = 0; j < h; ++j) {
    double i_sum = 0;
    for (int i = 0; i < w; ++i) {
      curr_pdf = pdf_envmap[j * w + i];
      marginal_distr += curr_pdf;
      if (prob_y[j] == 0)
        continue;
      i_sum += curr_pdf / prob_y[j];
      conds_y[(j + 1) * w + i] = i_sum;
    }
    marginal_y[j] = marginal_distr;
	}
	if (true)
		save_probability_debug();
}


// Helper functions

void EnvironmentLight::save_probability_debug() {
	uint32_t w = envMap->w, h = envMap->h;
	uint8_t* img = new uint8_t[4*w*h];

	for (int j = 0; j < h; ++j) {
		for (int i = 0; i < w; ++i) {
			img[4 * (j * w + i) + 3] = 255;
			img[4 * (j * w + i) + 0] = 255 * marginal_y[j];
			img[4 * (j * w + i) + 1] = 255 * conds_y[j * w + i];
		}
	}

    lodepng::encode("probability_debug.png", img, w, h);
    delete[] img;
}

Vector2D EnvironmentLight::theta_phi_to_xy(const Vector2D &theta_phi) const {
    uint32_t w = envMap->w, h = envMap->h;
    double x = theta_phi.y / 2. / M_PI * w;
    double y = theta_phi.x / M_PI * h;
    return Vector2D(x, y);
}

Vector2D EnvironmentLight::xy_to_theta_phi(const Vector2D &xy) const {
    uint32_t w = envMap->w, h = envMap->h;
    double x = xy.x;
    double y = xy.y;
    double phi = x / w * 2.0 * M_PI;
    double theta = y / h * M_PI;
    return Vector2D(theta, phi);
}

Vector2D EnvironmentLight::dir_to_theta_phi(const Vector3D &dir) const {
    dir.unit();
    double theta = acos(dir.y);
    double phi = atan2(-dir.z, dir.x) + M_PI;
    return Vector2D(theta, phi);
}

Vector3D EnvironmentLight::theta_phi_to_dir(const Vector2D& theta_phi) const {
    double theta = theta_phi.x;
    double phi = theta_phi.y;

    double y = cos(theta);
    double x = cos(phi - M_PI) * sin(theta);
    double z = -sin(phi - M_PI) * sin(theta);

    return Vector3D(x, y, z);
}

Spectrum EnvironmentLight::bilerp(const Vector2D& xy) const {
	uint32_t w = envMap->w;
	const std::vector<Spectrum>& data = envMap->data;
	double x = xy.x, y = xy.y;
	Spectrum ret;
	for (int i = 0; i < 4; ++i)
		ret += (i%2 ? x-floor(x) : ceil(x)-x) *
			   (i/2 ? y-floor(y) : ceil(y)-y) *
			   data[w * (floor(y) + i/2) + floor(x) + i%2];
	return ret;
}


Spectrum EnvironmentLight::sample_L(const Vector3D& p, Vector3D* wi,
                                    float* distToLight,
                                    float* pdf) const {
  // TODO: 3-2 Part 3 Tasks 2 and 3 (step 4)
	// First implement uniform sphere sampling for the environment light
	// Later implement full importance sampling

  // Part 3 Task 2 - uniform sampling
  // Vector3D sample_dir = sampler_uniform_sphere.get_sample();
  // Vector2D sample_xy = (dir_to_theta_phi(sample_dir));
  // Spectrum s = bilerp(sample_xy);
  // *distToLight = INF_D;
  // *pdf = 1.0 / (4.0 * PI);
  // return s;

  // Task 3, Step 4 - importance sampling
  uint32_t w = envMap->w, h = envMap->h;
  Vector2D sample2D = sampler_uniform2d.get_sample();
  std::vector<double>::iterator up_x, up_y;
  std::vector<double> row(marginal_y, marginal_y + h);
  up_y = std::upper_bound(row.begin(), row.end(), sample2D.y);
  int row_index = up_y - row.begin();

  std::vector<double> col(conds_y + w * row_index, conds_y + w * row_index + w - 1);
  up_x = std::upper_bound(col.begin(), col.end(), sample2D.x);
  int col_index = up_x - col.begin();

  Vector2D sample_xy = Vector2D(col_index, row_index);
  *wi = theta_phi_to_dir(xy_to_theta_phi(sample_xy));
  *distToLight = INF_D;
  // *pdf = pdf_envmap[int(sample_xy.x) + w * int(sample_xy.y)] * (w * h / (2 * PI * PI * sin(PI * int(sample_xy.y) / h)));
  *pdf = pdf_envmap[col_index + w * row_index] * (w * h / (2 * PI * PI * sin(PI * sample_xy.y / h)));
  // return bilerp(sample_xy);
  // return envMap->data[w * int(sample_xy.y) + int(sample_xy.x)];
  return envMap->data[w * row_index + col_index];
}

Spectrum EnvironmentLight::sample_dir(const Ray& r) const {
  // TODO: 3-2 Part 3 Task 1
	// Use the helper functions to convert r.d into (x,y)
	// then bilerp the return value
  Vector2D xy = theta_phi_to_xy(dir_to_theta_phi(r.d));
	return bilerp(xy);

}

} // namespace StaticScene
} // namespace CGL
