#pragma once
#include "integrator.h"

class UvIntegrator : public Integrator {
public:
  void render(const Scene& scene, Camera& cam, Image& out, IntegratorSettings& settings) {
    int xres = out.get_width();
    int yres = out.get_height();
    Intersection isect;

    for (int x = 0; x < xres; x++) {
      for (int y = 0; y < yres; y++) {
        Ray r = cam.getRay(x, y, xres, yres, false);

        bool hit = scene.intersect(r, 1.0, 9999999.0, isect);

        Color c = glm::vec3(0.0, 0.0, 0.0);
        if (hit) c = glm::vec3(isect.u);

        out.set_pixel(x, y, c);
      }
    }
  }
};