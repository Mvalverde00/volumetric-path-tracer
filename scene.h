#pragma once

#include <vector>

#include "geometry/geometry.h"
#include "medium/medium.h"
#include "bxdf/brdf.h"
#include "acceleration/acceleration.h"
#include "light/light.h"
#include "light/background_light.h"


class Scene {
  std::vector<Geometry*> objs;
  std::unordered_map<std::string, Brdf*> mats;
  std::vector<Light*> lights;
  BackgroundLight background_light;

  // TODO: For now only one medium, in the future we should support many
  Medium* medium;

  Acceleration* accel;

public:
  Scene();
  ~Scene();

  // When adding an object/mat/medium, you give full ownership of the data
  // to the Scene class.
  void add_geometry(Geometry* obj);
  void add_mat(std::string name, Brdf* mat);
  void set_medium(Medium* med);
  void set_background(BackgroundLight light);

  Brdf* get_mat(std::string name);

  void init_acceleration(Acceleration* acceleration);

  bool intersect(const Ray& r, float t_min, float t_max, Intersection& isect) const;

  bool intersect_media(const Ray &r, float t_min, float t_max, Intersection &isect) const;

  Color sample_lights(const Intersection &isect, const glm::vec3 wo) const;

  Color background_color(const glm::vec3& dir) const;

  Color transmission(const Ray &r, float t_max) const;
};
