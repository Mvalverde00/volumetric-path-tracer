#pragma once

#include "../geometry/geometry.h"
#include "../light/light.h"

class ShapeLight : public Light {
  Geometry *geom;
  Color color;

public:
  Color emitted(const glm::vec3 dir) const { return color; };
  void sample(glm::vec3 &p, glm::vec3 &n) { geom->sample(p, n); };
  float pdf(const glm::vec3 &lp, const glm::vec3 &p) {
    float dist_sqd = glm::dot(lp - p, lp - p);
    return geom->area() / dist_sqd;
  };

  ShapeLight(Geometry *g, Color c) : geom(g), color(c){};
};
