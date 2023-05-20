#include "scene.h"

#include <iostream>

Scene::Scene() : objs(), mats(), accel(NULL) {}


Scene::~Scene() {
  for (Geometry* obj : objs) {
    delete obj;
  }

  for (auto& pair : mats) {
    delete pair.second;
  }

  if (accel != NULL) delete accel;
}

// When adding an object or mat, you give full ownership of the data
// to the Scene class.
void Scene::add_geometry(Geometry* obj) {
  objs.push_back(obj);
}

void Scene::add_mat(std::string name, Brdf* mat) {
  if (mats.find(name) != mats.end()) {
    std::cout << "Two materials with same name '" << name << "'.\n";
    exit(1);
  }

  mats[name] = mat;
}

void Scene::set_medium(Medium* med) {
  medium = med;
}

Brdf* Scene::get_mat(std::string name) {
  if (mats.find(name) == mats.end()) {
    std::cout << "No material with name '" << name << "'.\n";
    exit(1);
  }

  return mats[name];
}

void Scene::init_acceleration(Acceleration* acceleration) {
  std::cout << "Scene has " << objs.size() << " objects and " << mats.size() << " materials\n";

  accel = acceleration;
  accel->build(objs);
}

bool Scene::intersect(const Ray& r, float t_min, float t_max, Intersection& isect) const {
  return accel->intersect(r, t_min, t_max, isect);
}

bool Scene::intersect_media(const Ray &r, float t_min, float t_max,
                      Intersection &isect) const {
  return medium->intersect(r, t_min, t_max, isect);
}
