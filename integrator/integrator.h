#pragma once

#include "../scene.h"
#include "../image.h"
#include "../camera.h"

struct IntegratorSettings {
  int max_depth;
  int samples;
  int threads;

  IntegratorSettings() : max_depth(5), samples(1), threads(1) {};

  IntegratorSettings(int max_depth, int samples, int threads)
      : max_depth(max_depth), samples(samples), threads(threads) {};
};

class Integrator {
public:
  virtual void render(const Scene& scene, Camera& cam, Image& out, IntegratorSettings& settings) = 0;

  //void render_tile(ImageTile tile, const Scene& scene, Camera& cam, Image& out, IntegratorSettings& settings);
};
