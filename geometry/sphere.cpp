#include "sphere.h"

#include "../acceleration/aabb.h"

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