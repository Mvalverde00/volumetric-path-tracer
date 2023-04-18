#pragma once

#include "../geometry/geometry.h"

/* Defines a set of structures which can be used to accelerate
 * ray intersections with geometry */
class Acceleration {
public:
  virtual void build(std::vector<Geometry*>& objs) = 0;
  virtual bool intersect(const Ray& r, float t_min, float t_max, Intersection& isect) = 0;
};
