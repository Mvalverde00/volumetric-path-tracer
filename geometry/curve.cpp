#include "curve.h"

#include "../acceleration/aabb.h"
#include "../util.h"

#include <glm/gtc/matrix_transform.hpp>
#include <algorithm>
#include <iostream>

const float inv_log_4 = 1.0f / std::log(4.0f);
const float PI = 3.14159274101257324219;

int Curve::max_depth() {
  int n = 3;

  float L0 = 0;
  for (int i = 0; i <= n-2; i++) {
    L0 = std::max(L0, (float)fabs(cp[i].x - 2.0f * cp[i + 1].x + cp[i + 2].x));
    L0 = std::max(L0, (float)fabs(cp[i].y - 2.0f * cp[i + 1].y + cp[i + 2].y));
    L0 = std::max(L0, (float)fabs(cp[i].z - 2.0f * cp[i + 1].z + cp[i + 2].z));
  }

  float eps = width / 20.0f;
  float r = (std::sqrt(2.0f) * n * (n - 1) * L0) / (8.0f * eps);

  return (int)(std::log(r) * inv_log_4);
}

void Curve::project(const Ray& r, const glm::mat3& proj_matrix, glm::vec3 transformed[4]) {
  for (int i = 0; i < 4; i++) {
    glm::vec3 p = cp[i] - r.origin;
    transformed[i] = glm::vec3(glm::dot(proj_matrix[0], p),
                               glm::dot(proj_matrix[1], p),
                               glm::dot(proj_matrix[2], p));
  }
}

Curve::Curve() : width(1.0), Geometry(NULL) {
  for (int i = 0; i < 4; i++) {
    cp[i] = glm::vec3(1.0f);
  }

  max_rec_depth = 2;
}

Curve::Curve(float width, const glm::vec3 c[4], Brdf* mat) : width(width), Geometry(mat) {
  for (int i = 0; i < 4; i++) {
    cp[i] = c[i];
  }

  max_rec_depth = std::max(0, std::min(max_depth(), 10));
};

inline glm::vec3 evaluate(const glm::vec3 cp[4], float t) {
  float sq = (1 - t) * (1 - t);
  return (1 - t) * sq * cp[0] + 3 * sq * t * cp[1]
         + 3 * (1 - t) * t * t * cp[2] + t * t * t * cp[3];
}

inline glm::vec3 tangent(const glm::vec3 cp[4], float t) {
  float sq = (1 - t) * (1 - t);
  return - 3*sq*cp[0] + 3*sq*cp[1] - 6*t*(1 - t)*cp[1]
         - 3*t*t*cp[2] + 6*t*(1 - t)*cp[2] + 3*t*t*cp[3];
}

AABB no_width_bbox(const glm::vec3 cp[4]) {
  glm::vec3 min = glm::min(glm::min(glm::min(cp[0], cp[1]), cp[2]), cp[3]);
  glm::vec3 max = glm::max(glm::max(glm::max(cp[0], cp[1]), cp[2]), cp[3]);

  return AABB(min, max);
}

inline glm::mat3 get_projection_matrix(const Ray& r) {
  const glm::vec3& l = r.dir;
  float d = std::sqrt(l.x * l.x + l.z * l.z);

  glm::vec3 lx, ly;
  if (d == 0.0f) {
    lx = glm::vec3(1.0, 0.0, 0.0);
    ly = glm::vec3(0.0, 0.0, -l.y);
  }
  else {
    lx = glm::vec3(l.z / d, 0.0, -l.x / d);
    ly = glm::vec3(lx.z * l.y, d, -l.y * lx.x);
  }

  return glm::mat3(lx, ly, l);
}

bool Curve::intersect(const Ray& r, float t_min, float t_max, Intersection& isect) {
  // Transform current curve by transformation matrix
  glm::mat3 proj_matrix = get_projection_matrix(r);
  glm::vec3 projected_c[4];
  project(r, proj_matrix, projected_c);

  float t = t_max;

  bool hit = converge(max_rec_depth, projected_c, 0, 1, t, isect);

  bool test = hit && (t_min < t) && (t < t_max);
  if (test) {
    isect.point = r.at(isect.t);

    glm::mat3 inv = glm::transpose(proj_matrix);

    isect.n = glm::vec3(glm::dot(inv[0], isect.n), glm::dot(inv[1], isect.n), glm::dot(inv[2], isect.n));
    isect.du = glm::vec3(glm::dot(inv[0], isect.du), glm::dot(inv[1], isect.du), glm::dot(inv[2], isect.du));
    isect.dv = glm::vec3(glm::dot(inv[0], isect.dv), glm::dot(inv[1], isect.dv), glm::dot(inv[2], isect.dv));
  }

  return hit && (t_min < t) && (t < t_max);
}

bool Curve::converge(int depth, glm::vec3 ccp[4], float v0, float vn, float& t, Intersection& isect) {
  float width1 = width / 2.0f;
  float width2 = width1 * width1;
  float eps = width / 20.0f;

  AABB b = no_width_bbox(ccp);

  // Test ray against bounding box.
  if (b.minimum.z >= t || b.maximum.z <= eps
    || b.minimum.x >= width1 || b.maximum.x <= -width1
    || b.minimum.y >= width1 || b.maximum.y <= -width1) {
    return false;
  }
  else if (depth == 0) {

    glm::vec2 dir = glm::vec2(ccp[3] - ccp[0]);

    glm::vec2 dp0 = glm::vec2(ccp[1] - ccp[0]);
    if (glm::dot(dir, dp0) < 0)
      dp0 *= -1;
    if (glm::dot(dp0, glm::vec2(-ccp[0])) < 0) {
      return false;
    }

    glm::vec2 dpn = glm::vec2(ccp[3] - ccp[2]);
    if (glm::dot(dir, dpn) < 0)
      dpn *= -1;
    if (glm::dot(dpn, glm::vec2(ccp[3])) < 0) {
      return false;
    }

    // Approximate with line
    float w = dir.x * dir.x + dir.y * dir.y;
    if (w == 0) 
      return false;

    w = -(ccp[0].x * dir.x + ccp[0].y * dir.y) / w;
    w = glm::clamp(w, 0.0f, 1.0f);

    // u runs along 1D bezier curve
    float u = v0 * (1.0f - w) + vn * w;

    glm::vec3 p = evaluate(ccp, w);
    float origin_dist_sq = p.x * p.x + p.y * p.y;
    if (origin_dist_sq > width2 || p.z <= eps) {
      return false;
    }
    if (t < p.z)
      return false;

    // Compute v coordinate which runs along width.
    glm::vec3 tan = tangent(ccp, w);
    float edgeFunc = tan.x * -p.y + p.x * tan.y;
    float origin_dist = std::sqrt(origin_dist_sq);
    //float v = (edgeFunc > 0) ? 0.5f - origin_dist / width1 :
    //                           0.5f + origin_dist / width1;

    /*
    float v = (edgeFunc > 0) ? 0.5f + origin_dist / width :
                               0.5f - origin_dist / width;
    
    // 90 degree rotation in the X-Y plane
    //glm::vec3 bitan = -glm::vec3(-tan.y, tan.x, 0);
    glm::vec3 bitan = glm::vec3(-tan.y, tan.x, 0);

    // Give cylindrical appearance by rotating bitan
    float rot = -PI / 2.0f + v * PI;
    glm::mat4 rotMat = glm::rotate(glm::mat4(1.0f), -rot, tan);
    bitan = glm::vec3(rotMat * glm::vec4(bitan, 0.0f));
    */
    float v = (edgeFunc > 0) ? 0.5f + origin_dist / width :
                               0.5f - origin_dist / width;

    // 90 degree rotation in the X-Y plane
    //glm::vec3 bitan = -glm::vec3(-tan.y, tan.x, 0);
    glm::vec3 bitan = glm::vec3(-tan.y, tan.x, 0);

    // Give cylindrical appearance by rotating bitan
    float rot = -PI / 2.0f + v * PI;
    glm::mat4 rotMat = glm::rotate(glm::mat4(1.0f), -rot, tan);
    bitan = glm::vec3(rotMat * glm::vec4(bitan, 0.0f));

    t = p.z;
    isect.t = t;
    isect.brdf = material;
    isect.u = u;
    isect.v = v;
    isect.n = glm::normalize(glm::cross(tan, bitan));
    isect.err = 2.0f * width;

    isect.du = glm::normalize(tan);
    isect.dv = glm::normalize(bitan);
    return true;
  }
  else {
    // Recursively split curve.

    depth--;
    float vm = (v0 + vn) / 2.0f;
    glm::vec3 c_left[4], c_right[4];
    split(ccp, c_left, c_right);

    return converge(depth, c_left, v0, vm, t, isect)
           || converge(depth, c_right, vm, vn, t, isect);
  }

}

const glm::vec3 ONES(1.0, 1.0, 1.0);

inline void Curve::bbox(AABB& box) const {
  glm::vec3 min = glm::min(glm::min(glm::min(cp[0], cp[1]), cp[2]), cp[3]);
  glm::vec3 max = glm::max(glm::max(glm::max(cp[0], cp[1]), cp[2]), cp[3]);

  const float w2 = width / 2.0f;
  box = AABB(min - ONES * w2,
             max + ONES * w2);
}
