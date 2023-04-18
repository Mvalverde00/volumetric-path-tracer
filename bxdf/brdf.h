#pragma once

#include <glm/glm.hpp>
#include <json.hpp>

struct Intersection;
class Scene;

using Color = glm::vec3;

class Brdf {
public:
  virtual Color eval(const glm::vec3& wo, const glm::vec3& wi, const Intersection& isect) = 0;
  virtual glm::vec3 sample(const Intersection& isect, const glm::vec3& wo) = 0;
  virtual float pdf(const glm::vec3& wo, const glm::vec3& wi, const Intersection& isect) = 0;

  // calculate eval/pdf explicitly, which may allow for simplifications
  virtual Color evalWithPdf(const glm::vec3& wo, const glm::vec3& wi, const Intersection& isect) = 0;
};
