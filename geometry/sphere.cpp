#include "sphere.h"

#include "../acceleration/aabb.h"
#include "../util.h"

bool Sphere::intersect(const Ray& r, float t_min, float t_max, Intersection& isect) {
  // Math taken from https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection

  glm::vec3 oc = r.origin - center;
  float a = glm::dot(r.dir, r.dir);
  float b = 2.0f * glm::dot(r.dir, oc);
  float c = glm::dot(oc, oc) - radius * radius;

  float discriminant = b * b - 4.0f * a * c;
  if (discriminant < 0.0f) {
    return false;
  }
  float sqrtd = std::sqrt(discriminant);

  float root = (-b - sqrtd) / (2.0f * a);
  if (root < t_min || root > t_max) {
    root = (-b + sqrtd) / (2.0f * a);
    if (root < t_min || root > t_max) {
      return false;
    }
  }

  isect.t = root;
  isect.point = r.at(root);
  isect.n = (isect.point - center) / radius;
  isect.brdf = material;
  isect.u = 0.0f;
  isect.v = 0.0f;
  isect.err = 0.001f;
  return true;
};

void Sphere::bbox(AABB& box) const {
  box.minimum = center - radius * glm::vec3(1.0, 1.0, 1.0);
  box.maximum = center + radius * glm::vec3(1.0, 1.0, 1.0);
}

void Sphere::sample(glm::vec3 & p, glm::vec3 & n) const {
  // Assume the [p] we are passed in initially represents the point on the surface
  // we are trying to cast a light ray from, and [n] is the corresponding normal.
  glm::vec3 isect_point = p;

  // Proceed normally
  n = randUniformSphere();
  p = radius * n + center;

  // Then we can always generate a point visible to the light, by sampling the
  // correct hemisphere.
  // NOTE: This means in [area], we only use the surface area of a hemisphere.
  glm::vec3 ldir = p - isect_point;
  if (glm::dot(ldir, n) > 0) {
    n = -n;
  }
  p = radius * n + center;
}

float Sphere::area() const {
  // Since we are only considering the relevant hemisphere in [sample],
  // this is the area of one hemisphere.
  return 2 * 3.14159265f * radius * radius;
}
