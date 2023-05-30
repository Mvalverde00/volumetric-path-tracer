#pragma once
#include "../ray.h"
#include "../bxdf/brdf.h"

#include <iostream>

class AABB;

struct Intersection {
  float t;
  glm::vec3 point;
  glm::vec3 n;

  // Necessary for texturing
  float u, v;

  // Magnitude of possible error in intersection point.
  // Should offset newly spawned ray by this distance to
  // avoid self-intersection.
  float err;

  // Derivatives of u and v in world coords.
  glm::vec3 du, dv;

  Brdf* brdf;
};

class Geometry {
public:
  Brdf *material;


  virtual bool intersect(const Ray& r, float t_min, float t_max, Intersection& isect) = 0;
  virtual void bbox(AABB& box) const = 0;

  virtual void sample(glm::vec3& p, glm::vec3& n) const {
    std::cout << "Sampling not implemented for this geometry\n";
    exit(1);
  };
  virtual float area() const {
    std::cout << "Area not implemented for this geometry\n";
    exit(1);
  };

  bool is_light() const {return material->emit() != Color(0.0, 0.0, 0.0);};

  Geometry(Brdf* mat) : material(mat) {};
};
