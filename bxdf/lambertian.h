#pragma once
#include "brdf.h"

#include "../util.h"
#include "../geometry/geometry.h"

const float PI = 3.141592f;

class Lambertian : public Brdf {
  Color diffuse;

public:
  Lambertian(Color c) : diffuse(c) {};

  Color eval(const glm::vec3& wo, const glm::vec3& wi, const Intersection& isect) {
    return diffuse * glm::dot(isect.n, wi) / PI;
  }

  std::optional<glm::vec3> sample(const Intersection& isect, const glm::vec3& wo) {
    // Sample wi uniformly from hemisphere aligned with
    // +z axis.
    //glm::vec3 wi = randUniformHemisphere();
    glm::vec3 wi = randCosineHemisphere();

    // Rotate wi to match actual normal orientation.
    return rotateVector(orientFromZAxis(isect.n), wi);
  }

  float pdf(const glm::vec3& wo, const glm::vec3& wi, const Intersection& isect) {
    //return 1.0 / (2.0 * PI); // Uniform Hemisphere
    return glm::dot(isect.n, wi) / PI; // Cosine Hemisphere
  }

  Color evalWithPdf(const glm::vec3& wo, const glm::vec3& wi, const Intersection& isect) {
    return diffuse * glm::dot(isect.n, wi) * 2.0f;
  }

  Lambertian(const nlohmann::json& desc) {
    diffuse = parse_vec3(desc.at("albedo"));
  }

};
