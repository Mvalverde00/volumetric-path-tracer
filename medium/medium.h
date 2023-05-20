#pragma once

#include "../geometry/geometry.h"
#include "../ray.h"

class Medium {
public:
  virtual bool intersect(const Ray &r, float t_min, float t_max,
                         Intersection &isect) = 0;
};
