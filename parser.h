#pragma once
#include <string>

/*
class Image;
class Camera;
class Scene;
struct IntegratorSettings;
class Integrator;
*/

#include "integrator/integrator.h"
#include "image.h"
#include "camera.h"
#include "scene.h"


// All the data necessary to represent and render
// a scene
struct RenderData {
  Image out;
  std::string out_file;

  Camera cam;
  Scene scene;
  IntegratorSettings isettings;
  Integrator* integrator;
};

// Parses the scene described by the given file and loads necessary
// data into the provided classes. Returns true if parsing/loading was
// successful and false if any errors were encountered. If errors were
// encountered, the provided classes may be uninitialized/half-initialized.
bool parse_scene(const std::string& fname, RenderData& isettings);