#include "parser.h"

#include "json.hpp"
#include "glm/glm.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <openvdb/openvdb.h>

#include "integrator/normal_integrator.h"
#include "integrator/naive_path_integrator.h"
#include "integrator/uv_integrator.h"

#include "acceleration/sah_bvh.h"
#include "acceleration/geometry_list.h"

#include "bxdf/lambertian.h"
#include "bxdf/marschner_hair.h"
#include "bxdf/light.h"

#include "geometry/sphere.h"
#include "geometry/curve.h"

#include "medium/grid.h"


using namespace nlohmann;

Acceleration* parse_acceleration(const std::string& acc_str) {
  if (acc_str == "none") {
    return new GeometryList();
  }
  else if (acc_str == "bvh") {
    return new SahBvh();
  }

  std::cout << "Error parsing acceleration: '" << acc_str << "' is not a valid acceleration structure.\n";
  exit(1);
}

Integrator* parse_integrator(const std::string& int_str) {
  if (int_str == "path") {
    return new NaivePathIntegrator();
  }
  else if (int_str == "normal") {
    return new NormalIntegrator();
  }
  else if (int_str == "uv") {
    return new UvIntegrator();
  }

  std::cout << "Error parsing integrator: '" << int_str << "' is not a valid integrator.\n";
  exit(1);
}

void parse_brdf(json desc, Scene& scene) {
  std::string type = desc.at("type");

  if (type == "lambert") {
    scene.add_mat(desc.at("name"), new Lambertian(desc));
    return;
  }
  else if (type == "hair") {
    scene.add_mat(desc.at("name"), new MarschnerHair(desc));
    return;
  }
  else if (type == "light") {
    scene.add_mat(desc.at("name"), new Light(desc));
    return;
  }

  std::cout << "Error parsing brdf: '" << type << "' is not a valid brdf.\n";
  exit(1);
}

void parse_media(json desc, Scene& scene) {
  std::string type = desc.at("type");

  if (type == "grid") {
    openvdb::initialize();

    openvdb::io::File file(desc.at("file"));

    file.open();
    scene.set_medium(new Grid(file));
    file.close();
    return;
  }

  std::cout << "Error parsing medium: '" << type << "' is not a valid medium .\n";
  exit(1);
}

void parse_hair_file(std::string filename, std::string matname, Scene& scene) {
  std::ifstream file(filename);

  if (!file.is_open()) {
    std::cout << "ERR: Could not load hair file '" << filename << "'\n";
    exit(1);
  }

  std::string dummy;
  while (file >> dummy) {

    // Read useless items at beginning;
    for (int i = 0; i < 9; i++) {
      file >> dummy;
    }

    // Read 4 control points in
    glm::vec3 control[4];
    for (int i = 0; i < 4; i++) {
      file >> control[i].x >> control[i].y >> control[i].z;
    }

    // Skip more filler text
    for (int i = 0; i < 4; i++) {
      file >> dummy;
    }

    // first width
    float width1, width2;
    file >> width1;

    // filler text
    for (int i = 0; i < 4; i++) {
      file >> dummy;
    }

    // second with
    file >> width2;

    // filler at end of line
    file >> dummy;

    // Their format supports variable width which is interpolated along curve.
    // Instead, we use uniform width.
    float width = (width1 + width2) / 2.0f;

    
    // Split curves into multiple pieces to get better
    // BVH performance.
    const int SPLITS = 3;
    const int num_curves = 1 << SPLITS;
    glm::vec3 final_control[num_curves][4];

    // Seed with starting curve
    for (int i = 0; i < 4; i++) {
      final_control[0][i] = control[i];
    }

    for (int i = 0; i < SPLITS; i++) {
      int start = (1 << i) - 1;
      for (int j = start; j >= 0; j--) {
        split(final_control[j], final_control[2*j], final_control[2*j + 1]);
      }
    }

    for (int i = 0; i < num_curves; i++) {
      scene.add_geometry(new Curve(width, final_control[i], scene.get_mat(matname)));
    }
    //scene.add_geometry(new Curve(width, control, scene.get_mat(matname)));
  }
}

void parse_obj(json desc, Scene& scene) {
  std::string prim = desc.at("type");

  if (prim == "sphere") {
    glm::vec3 c = parse_vec3(desc.at("center"));
    float r = desc.at("radius");
    std::string mat_name = desc.at("bsdf");
    scene.add_geometry(new Sphere(c, r, scene.get_mat(mat_name)));
  }
  else if (prim == "curves") {
    //scene.add_geometry(new Sphere(glm::vec3(-0.4336, 12.39, 0.1425), 1.0f, scene.get_mat("hair")));
    //scene.add_geometry(new Sphere(glm::vec3(0.1342, 9.92412, -0.828), 1.0f, scene.get_mat("hair")));
    parse_hair_file(desc.at("file"), desc.at("bsdf"), scene);
  }
  else {
    std::cout << "Error parsing primitive: '" << prim << "' is not a valid primitive.\n";
    exit(1);
  }
}


void parse_camera(json desc, RenderData& rd) {
  int xres = desc.at("xres");
  int yres = desc.at("yres");
  float fovy = desc.at("fov");

  json trans = desc.at("transform");

  glm::vec3 pos = parse_vec3(trans.at("position"));
  glm::vec3 target = parse_vec3(trans.at("look_at"));
  glm::vec3 up = parse_vec3(trans.at("up"));
  glm::mat4 rot = Camera::lookAt(pos, target, up);

  float aspect = float(xres) / yres;
  float near = 1.0f;
  Frustum f = Frustum(aspect, fovy, near);

  rd.cam = Camera(pos, rot, f);
  rd.out = Image(xres, yres);
  rd.out.init();
}


Acceleration* parse_renderer(json desc, RenderData& rd) {
  int bounces = desc.at("max_bounces");
  int spp = desc.at("spp");
  int threads = desc.at("threads");

  rd.isettings.max_depth = bounces;
  rd.isettings.samples = spp;
  rd.isettings.threads = threads;


  rd.out_file = desc.at("output_file");
  rd.integrator = parse_integrator(desc.at("type"));
  return parse_acceleration(desc.at("acceleration"));
}


bool parse_scene(const std::string& fname, RenderData& rd) {
  std::ifstream file(fname);

  if (!file.is_open()) {
    std::cout << "ERR: Could not load file " << fname << "\n";
    return false;
  }

  std::cout << "Parsing file\n";
  json js;
  try {
    js = json::parse(file);
  }
  catch (json::parse_error& ex) {
    std::cout << "Error parsing scene file:\n" << ex.what() << "\n";
    exit(1);
  }

  json camera = js.at("camera");
  json renderer = js.at("renderer");

  json mats = js.at("bsdfs");
  json objs = js.at("primitives");
  json media = js.at("media");

  // Parse render details and store acceleration structure for later usage
  std::cout << "Parsing renderer...\n";
  Acceleration* acc_struct = parse_renderer(renderer, rd);
  std::cout << "Parsing camera...\n";
  parse_camera(camera, rd);

  std::cout << "\nTarget image size: (" << rd.out.get_width() << ", "
            << rd.out.get_height() << ")\n";
  std::cout << "Rendering at " << rd.isettings.samples << " spp with "
            << rd.isettings.threads << " threads\n\n";

  std::cout << "Parsing materials...\n";
  // Load all bsdfs
  for (auto& mat : mats) {
    parse_brdf(mat, rd.scene);
  }

  std::cout << "Parsing objects...\n";
  // Finally, load scene objects and build acceleration
  for (auto& obj : objs) {
    parse_obj(obj, rd.scene);
  }

  for (auto &medium : media) {
    parse_media(medium, rd.scene);
  }

  std::cout << "Creating acceleration structure...\n";
  rd.scene.init_acceleration(acc_struct);

  return true;
}


