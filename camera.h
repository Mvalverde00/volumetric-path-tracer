#pragma once
#include <glm/glm.hpp>

#include "ray.h"

struct Frustum {
  float aspect;
  float fovy; // Stored in degrees
  float near;

  Frustum(float aspect, float fovy, float near);
  Frustum();
};

class Camera {
  glm::vec3 pos;
  glm::mat4 rot;
  Frustum frustum;

public:
  Camera(glm::vec3 pos, glm::mat4 rot, Frustum frustum);
  Camera(glm::vec3 pos, glm::mat4 rot);
  Camera(glm::vec3 pos);
  Camera();

  Ray getRay(int i, int j, int XRES, int YRES, bool random);

  static glm::mat4 lookAt(glm::vec3 eye, glm::vec3 center, glm::vec3 up);

};
