#pragma once

#include <glm/glm.hpp>

using Color = glm::vec3;
class Light {
public:
  virtual Color emitted(const glm::vec3 dir) const = 0;
  virtual void sample(glm::vec3 &p, glm::vec3 &n) = 0;
  virtual float pdf(const glm::vec3 &lp, const glm::vec3 &p) = 0;
};
