#pragma once
#include "brdf.h"

#include "../util.h"
#include "../geometry/geometry.h"

class MarschnerHair : public Brdf {
  Color sigma_a;
  float eta;

  float s;
  float v[4]; // Longitudinal Variance
  float sin_2k_alpha[3];
  float cos_2k_alpha[3];

  std::vector<Color> A_p(float h, float cos_theta_o, const Color& t);

  std::vector<float> get_ap_pdf(float h, float cos_theta_o);

public:
  Color eval(const glm::vec3& wo, const glm::vec3& wi, const Intersection& isect);

  std::optional<glm::vec3> sample(const Intersection& isect, const glm::vec3& wo);

  float pdf(const glm::vec3& wo, const glm::vec3& wi, const Intersection& isect);

  Color evalWithPdf(const glm::vec3& wo, const glm::vec3& wi, const Intersection& isect);

  MarschnerHair(const nlohmann::json& desc);
};
