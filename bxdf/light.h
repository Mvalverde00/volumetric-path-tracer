#pragma once
#include "brdf.h"

#include "../util.h"
#include "../geometry/geometry.h"

class LightBrdf : public Brdf {
  Color emitted;

public:
  LightBrdf(Color c) : emitted(c) {};

  Color eval(const glm::vec3& wo, const glm::vec3& wi, const Intersection& isect) {
    return Color(0.0, 0.0, 0.0);
  }

  std::optional<glm::vec3> sample(const Intersection& isect, const glm::vec3& wo) {
    return std::nullopt;
  }

  float pdf(const glm::vec3& wo, const glm::vec3& wi, const Intersection& isect) {
    return 0.0;
  }

  Color evalWithPdf(const glm::vec3& wo, const glm::vec3& wi, const Intersection& isect) {
    return Color(0.0, 0.0, 0.0);
  }

  Color emit() override {return emitted;}

  LightBrdf(const nlohmann::json& desc) {
    emitted = parse_vec3(desc.at("emitted"));
  }

};
