//#define GLM_FORCE_LEFT_HANDED

#include <chrono>
#include <iostream>

#include "acceleration/geometry_list.h"
#include "acceleration/sah_bvh.h"
#include "bxdf/lambertian.h"
#include "camera.h"
#include "geometry/curve.h"
#include "geometry/sphere.h"
#include "image.h"
#include "integrator/naive_path_integrator.h"
#include "integrator/normal_integrator.h"
#include "scene.h"

#include "parser.h"

#include "util.h"

#include "bxdf/marschner_hair.h"

float RANDOM_FLOAT() { return rand() / (float)RAND_MAX; }

int main(int argc, char *argv[]) {
  RenderData rd;

  auto start = std::chrono::system_clock::now();
  if (argc == 2) {
    std::string fname = std::string(argv[1]);
    if (!parse_scene(fname, rd)) {
      std::cout << "Encountered errors while parsing, quitting program.\n";
      return 1;
    }
  } else {
    std::cout << "Usage: render.exe <scene_filename>\n";
    return 1;
  }

  auto end = std::chrono::system_clock::now();
  std::cout
      << "parsing & acceleration runtime: "
      << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
      << "\n";

  std::cout << "Beginning render...\n";
  start = std::chrono::system_clock::now();
  rd.integrator->render(rd.scene, rd.cam, rd.out, rd.isettings);
  end = std::chrono::system_clock::now();

  std::cout
      << "time elapsed: "
      << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
      << "\n";
  rd.out.write_to_file(rd.out_file);
}
