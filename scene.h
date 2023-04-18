#pragma once

#include <vector>

#include "geometry/geometry.h"
#include "bxdf/brdf.h"
#include "acceleration/acceleration.h"


class Scene {
 
  std::vector<Geometry*> objs;
  std::unordered_map<std::string, Brdf*> mats;
  //std::vector<Brdf*> mats;

  Acceleration* accel;

public:
  Scene();
  ~Scene();

  // When adding an object or mat, you give full ownership of the data
  // to the Scene class.
  void add_geometry(Geometry* obj);
  void add_mat(std::string name, Brdf* mat);

  Brdf* get_mat(std::string name);

  void init_acceleration(Acceleration* acceleration);

  bool intersect(const Ray& r, float t_min, float t_max, Intersection& isect) const;
};