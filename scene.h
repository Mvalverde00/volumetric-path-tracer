#pragma once

#include <vector>

#include "geometry/geometry.h"
#include "medium/medium.h"
#include "bxdf/brdf.h"
#include "acceleration/acceleration.h"


class Scene {
 
  std::vector<Geometry*> objs;
  std::unordered_map<std::string, Brdf*> mats;

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

  Brdf* get_mat(std::string name);

  void init_acceleration(Acceleration* acceleration);

  bool intersect(const Ray& r, float t_min, float t_max, Intersection& isect) const;

  bool intersect_media(const Ray &r, float t_min, float t_max, Intersection &isect) const;
};
