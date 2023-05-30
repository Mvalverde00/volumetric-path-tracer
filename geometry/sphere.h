#pragma once
#include <glm/glm.hpp>

#include "geometry.h"

class Sphere : public Geometry {
  glm::vec3 center;
  float radius;

public:
  Sphere() : center(glm::vec3(0.0, 0.0, 0.0)), radius(1.0f), Geometry(NULL) {};
  Sphere(glm::vec3 center, float radius, Brdf* mat) : center(center), radius(radius), Geometry(mat) {};

  bool intersect(const Ray& r, float t_min, float t_max, Intersection& isect);

  void bbox(AABB& box) const;

  void sample(glm::vec3 &p, glm::vec3 &n) const override;
  float area() const override;
};
