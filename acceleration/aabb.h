#pragma once
#include <glm/glm.hpp>

#include <algorithm>

#include "../ray.h"
#include "../geometry/geometry.h"

class AABB {
public:
  glm::vec3 minimum;
  glm::vec3 maximum;

  AABB() : minimum(0, 0, 0), maximum(0, 0, 0) {};
  AABB(const glm::vec3& a, const glm::vec3& b) : minimum(a), maximum(b) {};
  AABB(AABB& a, AABB& b) : minimum(glm::min(a.minimum, b.minimum)), maximum(glm::max(a.maximum, b.maximum)) {};

  bool hit(const Ray& r, float t_min, float t_max) const {
    for (int dim = 0; dim < 3; dim++) {
      float invD = 1.0f / r.dir[dim];
      float t0 = (minimum[dim] - r.origin[dim]) * invD;
      float t1 = (maximum[dim] - r.origin[dim]) * invD;

      // swap
      if (invD < 0.0f) {
        float temp = t0;
        t0 = t1;
        t1 = temp;
      }

      t_min = t0 > t_min ? t0 : t_min;
      t_max = t1 < t_max ? t1 : t_max;
      if (t_max < t_min) {
        return false;
      }
    }
    return true;
  };

  bool hit_t(const Ray& r, float t_min, float t_max, float& t_hit) const {
    t_hit = t_min;
    for (int dim = 0; dim < 3; dim++) {
      float invD = 1.0f / r.dir[dim];
      float t0 = (minimum[dim] - r.origin[dim]) * invD;
      float t1 = (maximum[dim] - r.origin[dim]) * invD;

      // swap
      if (invD < 0.0f) {
        float temp = t0;
        t0 = t1;
        t1 = temp;
      }

      t_min = t0 > t_min ? t0 : t_min;
      t_max = t1 < t_max ? t1 : t_max;

      if (t_max < t_min) {
        return false;
      }

      t_hit = std::max(t_min, t_hit);
    }
    return true;
  };

  // surface area.  For a point-like bounding box
  // (i.e. with 0 volume), we define the surface area to be infinity.
  inline float area() {
    if (minimum == maximum) {
      //std::cout << "Creating infinite area\n";
      return 999999999999.0f; // Infinity
    }

    float s1 = maximum.x - minimum.x;
    float s2 = maximum.y - minimum.y;
    float s3 = maximum.z - minimum.z;
    return 2.f * (s1 * s2 + s2 * s3 + s3 * s1);
  }

  // Checks if the two AABB intersect on the given axis, and returns the result.
  // If an intersection on the axis does occur, the min/max coords of the
  // overlapping region are also returned.
  bool intersect_axis(const AABB& other, int axis, float& imin, float& imax) {
    if (minimum[axis] > other.maximum[axis] || other.minimum[axis] > maximum[axis]) {
      return false;
    }

    imin = std::max(minimum[axis], other.minimum[axis]);
    imax = std::min(maximum[axis], other.maximum[axis]);
    return true;
  }

  bool intersect(const AABB& other, AABB& result) {
    glm::vec3 new_min = glm::vec3(0, 0, 0);
    glm::vec3 new_max = glm::vec3(0, 0, 0);
    for (int axis = 0; axis < 3; axis++) {
      if (!intersect_axis(other, axis, new_min[axis], new_max[axis]))
        return false;
    }

    result.minimum = new_min;
    result.maximum = new_max;
    return true;
  }
};


inline bool box_compare(const Geometry& a, const Geometry& b, const int axis) {
  AABB box_a, box_b;
  a.bbox(box_a);
  b.bbox(box_b);

  return box_a.minimum[axis] < box_b.minimum[axis];
}

inline bool box_compare_x(const Geometry& a, const Geometry& b) { return box_compare(a, b, 0); }
inline bool box_compare_y(const Geometry& a, const Geometry& b) { return box_compare(a, b, 1); }
inline bool box_compare_z(const Geometry& a, const Geometry& b) { return box_compare(a, b, 2); }


inline bool box_compare_centroid(const Geometry* a, const Geometry* b, const int axis) {
  AABB box_a, box_b;
  a->bbox(box_a);
  b->bbox(box_b);

  // Can skip multiplying both sides by 0.5f since it does not affect inequality.
  return (box_a.minimum[axis] + box_a.maximum[axis]) < (box_b.minimum[axis] + box_b.maximum[axis]);
}

inline bool box_compare_centroid_x(const Geometry* a, const Geometry* b) { return box_compare_centroid(a, b, 0); }
inline bool box_compare_centroid_y(const Geometry* a, const Geometry* b) { return box_compare_centroid(a, b, 1); }
inline bool box_compare_centroid_z(const Geometry* a, const Geometry* b) { return box_compare_centroid(a, b, 2); }
