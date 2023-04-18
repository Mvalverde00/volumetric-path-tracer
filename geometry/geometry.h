#pragma once
#include "../ray.h"
#include "../bxdf/brdf.h"

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
  virtual bool intersect(const Ray& r, float t_min, float t_max, Intersection& isect) = 0;
  virtual void bbox(AABB& box) const = 0;
};