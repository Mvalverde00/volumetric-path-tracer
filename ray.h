#pragma once

#include <glm/glm.hpp>

struct Ray {
  glm::vec3 origin;
  glm::vec3 dir;


  Ray() : origin(glm::vec3(0.0, 0.0, 0.0)), dir(glm::vec3(0.0, 0.0, 1.0)) {}
  Ray(glm::vec3 origin, glm::vec3 dir) : origin(origin), dir(dir) {}
  glm::vec3 at(float t) const { return origin + t * dir; }
};