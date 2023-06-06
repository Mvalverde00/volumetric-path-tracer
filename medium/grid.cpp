#include "grid.h"
#include "../util.h"

#include <openvdb/Types.h>
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tools/Interpolation.h>

float Grid::grid_density(int x, int y, int z) {
  return sampler->isSample(openvdb::Coord(x, y, z));
}

float Grid::density(glm::vec3 p) {
  return sampler->wsSample(openvdb::Vec3d(p.x, p.y, p.z));
}

bool Grid::is_inside(glm::vec3 p) { return bbox.contains(p); }

bool Grid::intersect(const Ray &r, float t_min, float t_max,
                     Intersection &isect) {
  // If we miss the bounding box, we miss the grid
  float t = 0.0;
  if (!bbox.hit_t(r, t_min, t_max, t)) {
    return false;
  }

  // We have an intersection, and [t] is set to the beginning of the intersection
  // with the bounding box.
  do {
    t -= logf(1.0f - randFloat()) / max_density;
  } while (density(r.at(t)) <= randFloat() * max_density && t <= t_max &&
           is_inside(r.at(t)));

  if (t <= t_max && is_inside(r.at(t))) {
    isect.t = t;
    return true;
  }

  return false;
}

Color Grid::Tr(const Ray &r, float t_max) {
  // If we miss the bounding box, we miss the grid
  float t = 0.0;
  if (!bbox.hit_t(r, 0.0f, t_max, t)) {
    return Color(1.0, 1.0, 1.0);
  }

  // We have an intersection, and [t] is set to the beginning of the
  // intersection with the bounding box.
  Color tr = Color(1.0, 1.0, 1.0);
  do {
    t -= logf(1.0f - randFloat()) / max_density;
    float d = density(r.at(t));
    tr *= 1.0 - std::max(0.0f, d / max_density);
  } while (t <= t_max && is_inside(r.at(t)));

  float albedo = 0.5;
  return tr * albedo;
}

Grid::Grid(openvdb::io::File &file) {
  // Read the grid from file
  openvdb::GridBase::Ptr baseGrid;
  for (openvdb::io::File::NameIterator nameIter = file.beginName();
       nameIter != file.endName(); ++nameIter) {
    if (nameIter.gridName() == "density") {
      baseGrid = file.readGrid(nameIter.gridName());
    }
  }
  grid = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);

  // Scale appropriately
  // TODO: This should be read in somehow, eventually
  grid->setTransform(openvdb::math::Transform::createLinearTransform(1));

  // Setup a sampler to read densities at points
  sampler =
      std::unique_ptr<openvdb::tools::GridSampler<openvdb::FloatGrid,
                                                  openvdb::tools::BoxSampler>>(
          new openvdb::tools::GridSampler<openvdb::FloatGrid,
                                          openvdb::tools::BoxSampler>(*grid));

  // Setup bounding box for volume
  auto bbox = grid->evalActiveVoxelBoundingBox();
  openvdb::math::Vec3d min = grid->indexToWorld(bbox.min());
  openvdb::math::Vec3d max = grid->indexToWorld(bbox.max());
  this->bbox = AABB(glm::vec3(min.x(), min.y(), min.z()),
                    glm::vec3(max.x(), max.y(), max.z()));

  std::cout << "Bounding box:" << min << ", " << max << "\n";

  // Calculate max density
  openvdb::math::MinMax<float> min_max =
      openvdb::tools::minMax(grid->tree(), true);
  max_density = min_max.max();
  std::cout << "min/max density = " << min_max.min() << "/" << min_max.max() << "\n";

  // TODO: Actually load these from a file
  sigma_a = Color(0.5, 0.5, 0.5);
  sigma_s = Color(0.5, 0.5, 0.5);
}
