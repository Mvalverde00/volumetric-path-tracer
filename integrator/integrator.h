#pragma once

#include "../camera.h"
#include "../image.h"
#include "../scene.h"

struct IntegratorSettings {
  int max_depth;
  int samples;
  int threads;
  bool allow_light_sampling;

  IntegratorSettings()
      : max_depth(5), samples(1), threads(1), allow_light_sampling(true){};

  IntegratorSettings(int max_depth, int samples, int threads,
                     bool allow_light_sampling)
      : max_depth(max_depth), samples(samples), threads(threads),
        allow_light_sampling(allow_light_sampling){};
};

class Integrator {
public:
  virtual void render(const Scene &scene, Camera &cam, Image &out,
                      IntegratorSettings &settings) = 0;

  // void render_tile(ImageTile tile, const Scene& scene, Camera& cam, Image&
  // out, IntegratorSettings& settings);
};
