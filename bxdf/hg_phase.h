#pragma once
#include "brdf.h"

#include "../util.h"
#include "../geometry/geometry.h"

class HgPhase : public Brdf {
  float g;
  Color albedo;

  inline float hgphase(float cos_theta) {
    float denom = 1.0 + g * g + 2 * g * cos_theta;
    return (1.0 - g * g) / (4.0 * 3.141592653f * denom * sqrt(denom));
  }

public:
  HgPhase(float g, Color a) : g(g), albedo(a) {};

  Color eval(const glm::vec3& wo, const glm::vec3& wi, const Intersection& isect) {
    return albedo * hgphase(glm::dot(wi, wo));
  }

  std::optional<glm::vec3> sample(const Intersection& isect, const glm::vec3& wo) {
    return randUniformSphere();
  }

  float pdf(const glm::vec3& wo, const glm::vec3& wi, const Intersection& isect) {
    return 1.0 / (4.0 * 3.141592653f); // Uniform Sphere
  }

  Color evalWithPdf(const glm::vec3& wo, const glm::vec3& wi, const Intersection& isect) {
    return albedo * hgphase(glm::dot(wi, wo)) * 4.0f * 3.141592653f;
  }
};
