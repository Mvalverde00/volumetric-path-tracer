// Path Tracer V2.cpp : This file contains the 'main' function. Program execution begins and ends there.

//#define GLM_FORCE_LEFT_HANDED

#include <iostream>
#include <chrono>

#include "image.h"
#include "camera.h"
#include "scene.h"
#include "geometry/sphere.h"
#include "integrator/normal_integrator.h"
#include "integrator/naive_path_integrator.h"
#include "bxdf/lambertian.h"
#include "acceleration/geometry_list.h"
#include "acceleration/sah_bvh.h"
#include "geometry/curve.h"

#include "parser.h"

#include "util.h"

#include "bxdf/marschner_hair.h"

float RANDOM_FLOAT() { return rand() / (float)RAND_MAX; }

int main(int argc, char* argv[]) {
  RenderData rd;

  auto start = std::chrono::system_clock::now();
  if (argc == 2) {
    std::string fname = std::string(argv[1]);
    if (!parse_scene(fname, rd)) {
      std::cout << "Encountered errors while parsing, quitting program.\n";
      return 1;
    }
  }
  else {
    std::cout << "Usage: render.exe <scene_filename>\n";
    return 1;
  }

  auto end = std::chrono::system_clock::now();
  std::cout << "parsing & acceleration runtime: "
            << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
            << "\n";


  std::cout << "Beginning render...\n";
  start = std::chrono::system_clock::now();
  rd.integrator->render(rd.scene, rd.cam, rd.out, rd.isettings);
  end = std::chrono::system_clock::now();

  std::cout << "time elapsed: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << "\n";
  rd.out.write_to_file(rd.out_file);
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
