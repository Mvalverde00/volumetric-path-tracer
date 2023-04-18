#pragma once
#include "integrator.h"

#include <vector>
#include <thread>

#include "threadpool.h"

struct ImageTile {
  int x;
  int y;
  int width;
  int height;
};

class NaivePathIntegrator : public Integrator {

  Color ray_color(const Scene& scene, const Ray& r, int depth, float t_min) {
    // If maximum depth exceeded, pixel is black
    if (depth <= 0) return Color(0.0, 0.0, 0.0);

    Intersection isect;
    if (scene.intersect(r, t_min, 9999999.f, isect)) {
      glm::vec3 wo = -r.dir;
      glm::vec3 wi = isect.brdf->sample(isect, wo);
      float pdf = isect.brdf->pdf(wo, wi, isect);

      /*
      std::cout << "isect point = " << isect.point.x << ", " << isect.point.y << ", " << isect.point.z << "\n";
      std::cout << "isect v =" << isect.v << "\n";
      std::cout << "wo=" << wo.x << ", " << wo.y << ", " << wo.z << "\n";
      Color out = isect.brdf->eval(wo, wi, isect);
      std::cout << "pdf=" << pdf << ";  BRDF=" << out.x << "," << out.y << "," << out.z << "\n";
      std::cout << "wi=" << wi.x << ", " << wi.y << ", " << wi.z << "\n";
      */

      //return (isect.brdf->eval(wo, wi, isect) / pdf) * ray_color(scene, Ray(isect.point, wi), depth - 1, isect.err);
      
      // TODO: Doesn't work with Lambertian.
      return ((float)fabs(glm::dot(wi, isect.n))) * (isect.brdf->eval(wo, wi, isect) / pdf) * ray_color(scene, Ray(isect.point, wi), depth - 1, isect.err);
      
      // TODO: In theory this should be used as optimization, eventually.
      //return isect.brdf->evalWithPdf(wo, wi, isect) * ray_color(scene, Ray(isect.point, wi), depth - 1);
    }

    // Ray hits scene background
    return Color(0.7, 0.7, 1.0);
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

          // TODO: DELETE DEBUG
          //r.dir = glm::vec3(0.0116387, -0.0116387, -0.999865);

          // Reposition ray so that it begins on film plane
          r = Ray(r.at(1.0f), glm::normalize(r.dir));


          c += ray_color(scene, r, settings.max_depth, 0.001f);
          //std::cout << "===========\n";
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

    /*
    int xres = out.get_width();
    int yres = out.get_height();
    

    Intersection isect;
    for (int x = 0; x < xres; x++) {
      for (int y = 0; y < yres; y++) {
        Color c = glm::vec3(0.0, 0.0, 0.0);
        
        for (int sample = 0; sample < settings.samples; sample++) {
          Ray r = cam.getRay(x, y, xres, yres, true);
          // Reposition ray so that it begins on film plane
          r = Ray(r.at(1.0f), glm::normalize(r.dir));
          
          c += ray_color(scene, r, settings.max_depth);
        }

        out.set_pixel(x, y, c/float(settings.samples));
      }
    }*/
  }
};
