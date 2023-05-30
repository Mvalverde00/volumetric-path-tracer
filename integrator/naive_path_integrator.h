#pragma once
#include "integrator.h"
#include "../util.h"

#include <optional>
#include <vector>
#include <thread>

#include "threadpool.h"

struct ImageTile {
  int x;
  int y;
  int width;
  int height;
};

inline float hgphase(float g, float cos_theta) {
  float denom = 1.0 + g*g + 2*g*cos_theta;
  return (1.0 - g*g) / (4.0 * 3.141592653f * denom * sqrt(denom));
}

class NaivePathIntegrator : public Integrator {

  Color ray_color(const Scene& scene, const Ray& r, int depth, int max_depth, float t_min) {
    // If maximum depth exceeded, pixel is black
    if (depth <= 0) return Color(0.0, 0.0, 0.0);

    // Find intersection with surface, ignoring any media.
    Intersection isect, media_isect;
    bool has_intersection = scene.intersect(r, t_min, 9999999.f, isect);

    // Check where we would intersect media.
    bool media_intersection = scene.intersect_media(r, t_min, 9999999.f, media_isect);
    if (media_intersection && media_isect.t < isect.t) {
      float albedo = 0.5;
      // Isotropic phase function.
      glm::vec3 wo = -r.dir;
      glm::vec3 wi = randUniformSphere();
      float pdf = 1.0 / (4.0 * 3.141592653f);
      float phase = hgphase(-0.6, glm::dot(wi, wo));

      Color direct = Color(10.0, 1.0, 1.0); // Sample lights to estimate direct lighting.
      return direct + (albedo * phase / pdf) * ray_color(scene, Ray(r.at(media_isect.t), wi), depth - 1, max_depth, 0.0f);
    }
    else if (has_intersection)  {
      glm::vec3 wo = -r.dir;
      std::optional<glm::vec3> wi = isect.brdf->sample(isect, wo);

      if (wi.has_value()) {
        float pdf = isect.brdf->pdf(wo, *wi, isect);
        Color direct = scene.sample_lights(isect, wo);
        return direct + (isect.brdf->eval(wo, *wi, isect) / pdf) *
                         ray_color(scene, Ray(isect.point, *wi), depth - 1,
                                   max_depth, isect.err);
      } else {
        // glm::vec3 emitted = isect.brdf->emit();
        // return isect.brdf->emit();

        if (depth == max_depth) {
        glm::vec3 emitted = isect.brdf->emit();
        return isect.brdf->emit();
        } else {
          return Color(0.0, 0.0, 0.0);
        }
      }
    }

    // Ray hits scene background
    // return Color(0.7, 0.7, 1.0);
    return Color(0.0, 0.0, 0.0);
  }

  void render_tile(ImageTile tile, const Scene& scene, Camera& cam, Image& out, IntegratorSettings& settings) {
    int xres = out.get_width();
    int yres = out.get_height();

    for (int x = 0; x < tile.width; x++) {
      for (int y = 0; y < tile.height; y++) {
        int screenX = tile.width * tile.x + x;
        int screenY = tile.height * tile.y + y;

        if (screenX >= xres || screenY >= yres) continue;

        Color c = glm::vec3(0.0, 0.0, 0.0);
        for (int sample = 0; sample < settings.samples; sample++) {
          Ray r = cam.getRay(screenX, screenY, xres, yres, true);

          // Reposition ray so that it begins on film plane
          r = Ray(r.at(1.0f), glm::normalize(r.dir));


          c += ray_color(scene, r, settings.max_depth, settings.max_depth, 0.001f);
        }
        out.set_pixel(screenX, screenY, c / float(settings.samples));
      }
    }
  }


public:
  void render(const Scene& scene, Camera& cam, Image& out, IntegratorSettings& settings) {
    ThreadPool pool(settings.threads);

    int tileWidth = 16;
    int tileHeight = 16;

    for (int tileX = 0; tileX < ((out.get_width() - 1) / tileWidth) + 1; tileX++) {
      for (int tileY = 0; tileY < ((out.get_height() - 1) / tileHeight) + 1; tileY++) {
        ImageTile tile = { tileX, tileY, tileWidth, tileHeight };

        pool.enqueue([=, &scene, &cam, &out, &settings]() {
          render_tile(tile, scene, cam, out, settings);
        });
      }
    }
  }
};
