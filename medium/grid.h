#include "../acceleration/aabb.h"
#include "medium.h"
#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>

class Grid : public Medium {
  openvdb::FloatGrid::Ptr grid;
  std::unique_ptr<openvdb::tools::GridSampler<openvdb::FloatGrid,
                                              openvdb::tools::BoxSampler>>
      sampler;
  Color sigma_a, sigma_s;
  float max_density;
  AABB bbox;

  // Calculate density at grid cell
  float grid_density(int x, int y, int z);

  // Calculate density at any point, interpolating grid cell densities
  float density(glm::vec3 p);

  // Whether a point is inside the grid
  bool is_inside(glm::vec3 p);

public:
  Grid(openvdb::io::File &file);

  bool intersect(const Ray &r, float t_min, float t_max, Intersection &isect);
  Color Tr(const Ray &r, float t_max);
};
