#pragma once

#include "acceleration.h"

class GeometryList : public Acceleration {
  std::vector<Geometry*>* objs;

public:
  void build(std::vector<Geometry*>& objects) { objs = &objects; };
  bool intersect(const Ray& r, float t_min, float t_max, Intersection& isect) {
    float closest_t = t_max + 1.0f;
    bool scene_hit = false;
    Intersection closest;
    for (Geometry* obj : *objs) {
      bool hit = obj->intersect(r, t_min, t_max, isect);

      if (hit && isect.t < closest_t) {
        closest = isect;
        closest_t = isect.t;
        scene_hit = true;
      }
    }

    if (!scene_hit) return false;

    isect = closest;
    return true;
  }
};