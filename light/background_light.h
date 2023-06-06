#pragma once

#include "../geometry/geometry.h"
#include "../light/light.h"

#include "../util.h"

class BackgroundLight : public Light {
  Color color;

public:
  Color emitted(const glm::vec3 dir) const { return color; };
  void sample(glm::vec3 &p, glm::vec3 &n) {
    glm::vec3 dir = randUniformSphere();
    p = dir * 9999999.0f;
    n = -dir;
  };
  float pdf(const glm::vec3& lp, const glm::vec3& p) { return 1.0 / (4.0f * 3.14159265f);};

  BackgroundLight(Color c) : color(c){};
  BackgroundLight() : color(Color(0.0, 0.0, 0.0)) {};
};
