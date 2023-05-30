#pragma once

#include "../geometry/geometry.h"
#include "../light/light.h"


class ShapeLight : public Light {
  Geometry* geom;
  Color color;

public:
  Color emitted(const glm::vec3 dir) {return color;};
  void sample(glm::vec3 &p, glm::vec3 &n) {geom->sample(p, n);};
  float area() {return geom->area();};

  ShapeLight(Geometry* g, Color c) : geom(g), color(c) {};
};
