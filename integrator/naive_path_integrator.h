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

  Color ray_color(const Scene& scene, const Ray& r, int depth, float t_min) {
    // If maximum depth exceeded, pixel is black
    if (depth <= 0) return Color(0.0, 0.0, 0.0);

    // Find intersection with surface, ignoring any media.
    Intersection isect;
    bool has_intersection = scene.intersect(r, t_min, 9999999.f, isect);

    // Check where we would intersect media.
    float extinction_coeff = 0.4;
    float collision_t = -log(1.0f - randFloat()) / extinction_coeff;
    float albedo = 0.5;

    if (collision_t < isect.t) {
      // Isotropic phase function.
      glm::vec3 wo = -r.dir;
      glm::vec3 wi = randUniformSphere();
      float pdf = 1.0 / (4.0 * 3.141592653f);
      float phase = hgphase(-0.6, glm::dot(wi, wo));

      return (albedo * phase / pdf) * ray_color(scene, Ray(r.at(collision_t), wi), depth - 1, 0.0f);
    }
    else if (has_intersection)  {
      glm::vec3 wo = -r.dir;
      std::optional<glm::vec3> wi = isect.brdf->sample(isect, wo);

      if (wi.has_value()) {
        float pdf = isect.brdf->pdf(wo, *wi, isect);
        return (isect.brdf->eval(wo, *wi, isect) / pdf) *
               ray_color(scene, Ray(isect.point, *wi), depth - 1, isect.err);

        // TODO: Doesn't work with Lambertian.
        // return ((float)fabs(glm::dot(wi, isect.n))) * (isect.brdf->eval(wo,
        // wi, isect) / pdf) * ray_color(scene, Ray(isect.point, wi), depth - 1,
        // isect.err);

        // TODO: In theory this should be used as optimization, eventually.
        // return isect.brdf->evalWithPdf(wo, wi, isect) * ray_color(scene,
        // Ray(isect.point, wi), depth - 1);
      } else {
        glm::vec3 emitted = isect.brdf->emit();
        return isect.brdf->emit();
      }
    }

    // Ray hits scene background
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


          c += ray_color(scene, r, settings.max_depth, 0.001f);
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
