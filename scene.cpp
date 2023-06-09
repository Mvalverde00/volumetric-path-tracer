#include "scene.h"

#include "light/shape_light.h"

#include <iostream>

Scene::Scene() : objs(), mats(), accel(NULL), medium(NULL), lights() {}


Scene::~Scene() {
  for (Geometry* obj : objs) {
    delete obj;
  }

  for (auto& pair : mats) {
    delete pair.second;
  }

  if (accel != NULL) delete accel;

  if (medium != NULL) delete medium;

}

// When adding an object or mat, you give full ownership of the data
// to the Scene class.
void Scene::add_geometry(Geometry* obj) {
  objs.push_back(obj);

  if (obj->is_light()) {
    lights.push_back(new ShapeLight(obj, obj->material->emit()));
  }
}

void Scene::add_mat(std::string name, Brdf* mat) {
  if (mats.find(name) != mats.end()) {
    std::cout << "Two materials with same name '" << name << "'.\n";
    exit(1);
  }

  mats[name] = mat;
}

void Scene::set_medium(Medium *med) { medium = med; }

void Scene::set_background(BackgroundLight light) {
  background_light = light;

  // TODO: Currently background lights are all uniform color, so we can
  // plug in any vector here. We just want to ignore a pure-black background,
  // since summing over it in the direct light sampling is a waste.
  // Once environment maps are added, this will need some revision
  if (background_light.emitted(glm::vec3(0.0, 0.0, 0.0)) != Color(0.0, 0.0, 0.0)) {
    std::cout << "Preparing background for explicit light sampling...\n";
    lights.push_back(&background_light);
  }
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
  if (medium == NULL) return false;

  return medium->intersect(r, t_min, t_max, isect);
}

Color Scene::sample_lights(const Intersection &isect, const glm::vec3 wo) const {
  Color accum = Color(0.0, 0.0, 0.0);
  Intersection light_isect;
  glm::vec3 lp, ln; // light point and normal
  for (const auto light : lights) {
    lp = isect.point;
    light->sample(lp, ln);
    glm::vec3 ldir = glm::normalize(lp - isect.point);
    Ray r(isect.point, ldir);
    float n_dot_l = glm::dot(isect.n, ldir);
    float ln_dot_l = glm::dot(ln, -ldir);

    float dist = glm::length(lp - isect.point);
    bool blocked = intersect(r, isect.err, dist, light_isect);
    if (n_dot_l > 0 && ln_dot_l > 0){
      if (blocked && light_isect.t < dist - 0.001) {continue;}
      float pdf = light->pdf(lp, isect.point);
      Color tr = transmission(Ray(isect.point, ldir), dist);
      accum += tr * pdf * ln_dot_l * isect.brdf->eval(wo, ldir, isect) * light->emitted(-ldir);
      continue;
    }
  }
  return accum;
}

Color Scene::background_color(const glm::vec3& dir) const {
  return background_light.emitted(dir);
}

Color Scene::transmission(const Ray &r, float t_max) const {
  // TODO: Eventually we will loop over all media and multiply transmissions
  // from all of them
  if (medium) {
    return medium->Tr(r, t_max);
  }
  return Color(1.0, 1.0, 1.0);
}
