#pragma once

#include <vector>

#include "acceleration.h"
#include "aabb.h"

/* Implemenation of Surface Area Heuristic BVH taken from
 * https://graphics.stanford.edu/~boulos/papers/togbvh.pdf */

struct BVHNode {
  AABB bbox; // 6 floats = 24 bytes

  /* when positive, left/right store the indices of the
   * corresponding left/right BVHnodes.  When negative,
   * the BVHNode is a leaf node and the interval
   * [-left, -right) are the geometries contained. */
  int left; // 4 bytes
  int right; // 4 bytes
};

namespace BVH {
  // Create a leaf node corresponding to triangle indices in [start, end)
  BVHNode makeLeaf(AABB& bbox, int start, int end);

  // Create an internal node with the corresponding left and right nodes
  BVHNode makeInternal(AABB& bbox, int left, int right);

  // Abstract away common actions since BVH will likely undergo significant changes
  bool isLeaf(std::vector<BVHNode>& nodes, int idx);
  int getLeft(std::vector<BVHNode>& nodes, int idx);
  int getRight(std::vector<BVHNode>& nodes, int idx);
  int getStart(std::vector<BVHNode>& nodes, int idx);
  int getEnd(std::vector<BVHNode>& nodes, int idx);
  AABB& getAABB(std::vector<BVHNode>& nodes, int idx);

};

struct GeomWithAABB {
  Geometry* geom;
  AABB aabb;
};

class SahBvh : public Acceleration {
  std::vector<BVHNode> nodes;
  std::vector<Geometry*>* objs;

  int buildHelper(std::vector<GeomWithAABB>& objs, int start, int end, int depth);

public:
  SahBvh() : nodes() {};

  void build(std::vector<Geometry*>& objects);
  bool intersect(const Ray& r, float t_min, float t_max, Intersection& isect);
};
