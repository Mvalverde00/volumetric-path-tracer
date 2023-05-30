#pragma once

#include "geometry.h"

class Curve : public Geometry {
  float width;
  int max_rec_depth;
  // Control points for curve
  glm::vec3 cp[4];

  int max_depth();
  void project(const Ray& r, const glm::mat3& proj_matrix, glm::vec3 transformed[4]);

  bool converge(int depth, glm::vec3 ccp[4], float v0, float vn, float& t, Intersection& isect);

public:
  Curve();
  Curve(float width, const glm::vec3 cp[4], Brdf* mat);

  bool intersect(const Ray& r, float t_min, float t_max, Intersection& isect);
  inline void bbox(AABB& box) const;
};



// Implementation from https://stackoverflow.com/questions/18655135/divide-bezier-curve-into-two-equal-halves/18681336#18681336
inline void split(const glm::vec3 cp[4], glm::vec3 cp1[4], glm::vec3 cp2[4]) {
  glm::vec3 E = (cp[0] + cp[1]) / 2.0f;
  glm::vec3 F = (cp[0] + 2.0f * cp[1] + cp[2]) / 4.0f;
  glm::vec3 G = (cp[0] + 3.0f * cp[1] + 3.0f * cp[2] + cp[3]) / 8.0f;

  glm::vec3 H = (cp[1] + 2.0f * cp[2] + cp[3]) / 4.0f;
  glm::vec3 J = (cp[2] + cp[3]) / 2.0f;

  glm::vec3 left[4] = { cp[0], E, F, G };
  glm::vec3 right[4] = { G, H, J, cp[3] };

  for (int i = 0; i < 4; i++) {
    cp1[i] = left[i];
    cp2[i] = right[i];
  }
}
