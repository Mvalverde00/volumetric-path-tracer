#include "sah_bvh.h"

#include <algorithm>
#include <iostream>

BVHNode BVH::makeLeaf(AABB& bbox, int start, int end) {
  // Negate indices to be able to tell it is a leaf.
  return BVHNode{ bbox, -start, -end };
}

BVHNode BVH::makeInternal(AABB& bbox, int left, int right) {
  return BVHNode{ bbox, left, right };
}

bool BVH::isLeaf(std::vector<BVHNode>& nodes, int idx) {
  return getLeft(nodes, idx) < 0 || getRight(nodes, idx) < 0;
}

int BVH::getLeft(std::vector<BVHNode>& nodes, int idx) {
  return nodes[idx].left;
}
int BVH::getRight(std::vector<BVHNode>& nodes, int idx) {
  return nodes[idx].right;
}

AABB& BVH::getAABB(std::vector<BVHNode>& nodes, int idx) {
  return nodes[idx].bbox;
};

int BVH::getStart(std::vector<BVHNode>& nodes, int idx) {
  return -nodes[idx].left;
}
int BVH::getEnd(std::vector<BVHNode>& nodes, int idx) {
  return -nodes[idx].right;
}




const float C_tri = 1.0f;
const float C_step = 0.0f;

inline float SAH_Heuristic(int n_s1, float area_s1, int n_s2, float area_s2, float area_total) {
  return (2.f * C_step) + (area_s1 / area_total * n_s1 * C_tri) + (area_s2 / area_total * n_s2 * C_tri);
}

inline float SAH_Heuristic2(int n_s1, float area_s1, int n_s2, float area_s2) {
  // Since 2.0f * C_step is a constant, we can simply minimize the rest of the equation.
  return area_s1 * n_s1 + area_s2 * n_s2;
}

/* Optimized versions to run with precomputed AABBs*/
inline bool box_compare_centroid_opt(const GeomWithAABB& a, const GeomWithAABB& b, const int axis) {

  // Can skip multiplying both sides by 0.5f since it does not affect inequality.
  return (a.aabb.minimum[axis] + a.aabb.maximum[axis]) < (b.aabb.minimum[axis] + b.aabb.maximum[axis]);
}
inline bool box_compare_centroid_x_opt(const GeomWithAABB& a, const GeomWithAABB& b) { return box_compare_centroid_opt(a, b, 0); }
inline bool box_compare_centroid_y_opt(const GeomWithAABB& a, const GeomWithAABB& b) { return box_compare_centroid_opt(a, b, 1); }
inline bool box_compare_centroid_z_opt(const GeomWithAABB& a, const GeomWithAABB& b) { return box_compare_centroid_opt(a, b, 2); }

int SahBvh::buildHelper(std::vector<GeomWithAABB>& geom_aabb, int start, int end, int depth) {
  int span = end - start;
  float bestCost = C_tri * span;
  int bestAxis = -1;
  int bestEvent = -1;

  for (int axis = 0; axis < 3; axis++) {
    /*
    auto comparator = axis == 0 ? box_compare_centroid_x
      : axis == 1 ? box_compare_centroid_y
      : box_compare_centroid_z;
      */
    auto comparator = axis == 0 ? box_compare_centroid_x_opt
      : axis == 1 ? box_compare_centroid_y_opt
      : box_compare_centroid_z_opt;

    //std::sort(objs.begin() + start, objs.begin() + end, comparator);
    //std::cout << start << ", " << end << ", " << geom_aabb.size() << "\n";
    std::sort(geom_aabb.begin() + start, geom_aabb.begin() + end, comparator);

    std::vector<float> leftAreas(span);
    AABB totalBBox, objBBox;
    totalBBox = AABB();
    for (int i = start; i < end; i++) {
      leftAreas[i - start] = totalBBox.area();

      // Expand bbox to include tri
      //objs[i]->bbox(objBBox);
      objBBox = geom_aabb[i].aabb;
      if (i == start)
        totalBBox = objBBox;
      else
        totalBBox = AABB(totalBBox, objBBox);

    }
    float totalArea = totalBBox.area();

    std::vector<float> rightAreas(span);
    totalBBox = AABB(); // reset overall bbox
    for (int i = end - 1; i >= start; i--) {
      rightAreas[i - start] = totalBBox.area();

      // Expand bbox to include tri
      //objs[i]->bbox(objBBox);
      objBBox = geom_aabb[i].aabb;
      if (i == end - 1)
        totalBBox = objBBox;
      else
        totalBBox = AABB(totalBBox, objBBox);

      float thisCost = SAH_Heuristic(i - start + 1, leftAreas[i - start], span - i + start - 1, rightAreas[i - start], totalArea);

      if (thisCost < bestCost) {
        bestCost = thisCost;
        bestEvent = i;
        bestAxis = axis;
      }
    }

  }

  if (bestAxis == -1) {
    // Best option is to just make current selection a leaf node
    AABB objBBox;
    //objs[start]->bbox(objBBox);
    objBBox = geom_aabb[start].aabb;
    AABB total = objBBox;
    for (int i = start; i < end; i++) {
      //objs[i]->bbox(objBBox);
      objBBox = geom_aabb[i].aabb;
      total = AABB(total, objBBox);
    }

    nodes.push_back(BVH::makeLeaf(total, start, end));
    return int(nodes.size()) - 1;
  }
  else {
    /*
    auto comparator = bestAxis == 0 ? box_compare_centroid_x
      : bestAxis == 1 ? box_compare_centroid_y
      : box_compare_centroid_z;
      */
    auto comparator = bestAxis == 0 ? box_compare_centroid_x_opt
      : bestAxis == 1 ? box_compare_centroid_y_opt
      : box_compare_centroid_z_opt;

    //std::sort(objs.begin() + start, objs.begin() + end, comparator);
    std::sort(geom_aabb.begin() + start, geom_aabb.begin() + end, comparator);

    // reserve space for current internal node
    nodes.push_back(BVHNode());
    int idx = int(nodes.size()) - 1;

    //int left = buildHelper(objs, start, bestEvent, depth + 1);
    //int right = buildHelper(objs, bestEvent, end, depth + 1);
    int left = buildHelper(geom_aabb, start, bestEvent, depth + 1);
    int right = buildHelper(geom_aabb, bestEvent, end, depth + 1);

    // Update internal node with correct data
    AABB merged = AABB(BVH::getAABB(nodes, left), BVH::getAABB(nodes, right));
    nodes[idx] = BVH::makeInternal(merged, left, right);

    return idx;
  }
}

void SahBvh::build(std::vector<Geometry*>& objects) {
  std::vector<GeomWithAABB> geom_aabb(objects.size());
  for (int i = 0; i < objects.size(); i++) {
    geom_aabb[i].geom = objects[i];
    objects[i]->bbox(geom_aabb[i].aabb);
  }


  int _ = buildHelper(geom_aabb, 0, objects.size(), 0);


  for (int i = 0; i < objects.size(); i++) {
    objects[i] = geom_aabb[i].geom;
  }

  objs = &objects;
}

bool SahBvh::intersect(const Ray& r, float t_min, float t_max, Intersection& isect) {
  int stack[64];
  int* stackPtr = stack;
  *stackPtr++ = 0;

  Intersection temp;
  isect.t = t_max;
  bool hit = false;
  volatile int stub = 0;

  int node = 0;
  // A very small mesh might only have a leaf, so we must account
  // for that.
  if (BVH::isLeaf(nodes, node)) {
    for (int i = BVH::getStart(nodes, node); i < BVH::getEnd(nodes, node); i++) {
      if ((*objs)[i]->intersect(r, t_min, t_max, temp) && temp.t < isect.t) {
        isect = temp;
        hit = true;
      }
    }
    return hit;
  }

  do {
    int childL = BVH::getLeft(nodes, node);
    int childR = BVH::getRight(nodes, node);

    bool overlapL = BVH::getAABB(nodes, childL).hit(r, t_min, t_max);
    bool overlapR = BVH::getAABB(nodes, childR).hit(r, t_min, t_max);

    // Check for overlap against leaf nodes and add to list of
    // possible intersections if appropriate
    if (overlapL && BVH::isLeaf(nodes, childL)) {
      for (int i = BVH::getStart(nodes, childL); i < BVH::getEnd(nodes, childL); i++) {
        if ((*objs)[i]->intersect(r, t_min, t_max, temp) && temp.t < isect.t) {
          isect = temp;
          hit = true;
        }
      }
    }
    if (overlapR && BVH::isLeaf(nodes, childR)) {
      for (int i = BVH::getStart(nodes, childR); i < BVH::getEnd(nodes, childR); i++) {
        if ((*objs)[i]->intersect(r, t_min, t_max, temp) && temp.t < isect.t) {
          isect = temp;
          hit = true;
        }
      }
    }

    bool traverseL = overlapL && !BVH::isLeaf(nodes, childL);
    bool traverseR = overlapR && !BVH::isLeaf(nodes, childR);

    // If we have no active children left to traverse,
    // pop a node from the stack and start traversing it instead.
    if (!traverseL && !traverseR)
      node = *--stackPtr;
    else {
      // We have at least one active child to traverse.
      // Prioritize traversing the left child.  Push right
      // child to stack to traverse later, if applicable.
      node = traverseL ? childL : childR;
      if (traverseL && traverseR)
        *stackPtr++ = childR;
    }

  } while (node != 0);

  return hit || stub;
}
