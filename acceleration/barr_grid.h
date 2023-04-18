#pragma once

/* Implemenation of 3D Grid structure taken from
 * https://www.microsoft.com/en-us/research/wp-content/uploads/2017/01/p119-snyder.pdf
 */

/*
#include "aabb.h"
#include "acceleration.h"

const int BARR_OBJECT = 0;
const int BARR_LIST = 1;
const int BARR_GEOMETRY = 2;

struct BarrObject {
  AABB bbox;
  void* root_object;
  int object_type;

  BarrObject() {}
  BarrObject(AABB bbox, void* root, int type) : bbox(bbox), root_object(root), object_type(type) {}
};

struct BarrList {
  std::vector<BarrObject*> list;

};

class BarrGrid : public Acceleration {
  AABB bbox;
  int x_divs, y_divs, z_divs, n_cells;
  BarrObject** cells;

public:
  BarrGrid(AABB box, int x_divs, int y_divs, int z_divs) : bbox(box), x_divs(x_divs), y_divs(y_divs), z_divs(z_divs){
    n_cells = x_divs * y_divs * z_divs;
    cells = new BarrObject*[n_cells];


    glm::vec3 cell_step = (bbox.maximum - bbox.minimum);
    cell_step.x /= x_divs;
    cell_step.y /= y_divs;
    cell_step.z /= z_divs;
    for (int x = 0; x < x_divs; x++) {
      for (int y = 0; y < y_divs; y++) {
        for (int z = 0; z < z_divs; z++) {
          int idx = triple_to_idx(x, y, z);
          cells[idx] = new BarrObject();

          glm::vec3 minimum = bbox.minimum + cell_step * glm::vec3(x, y, z);
          glm::vec3 maximum = minimum + cell_step;
          cells[idx]->bbox = AABB(minimum, maximum);
        }
      }
    }
  }

  ~BarrGrid() {
    for (int i = 0; i < n_cells; i++) {
      // TODO: first, delete any objects inside the cells
      delete cells[i];
    }

    delete[] cells;
  }

  int triple_to_idx(int x, int y, int z) {
    return x + y * x_cells + z * x_cells * y_cells;
  }


  void build(std::vector<Geometry*>& objs) {
    for (Geometry* obj : objs) {
      AABB obj_box, intersection;
      obj->bbox(obj_box);
      for (int i = 0; i < x_divs * y_divs * z_divs; i++) {
        BarrObject* cell = cells[i];

        if (obj_box.intersect(cell->bbox, intersection)) {
          BarrObject* new_obj = new BarrObject(intersection, obj, BARR_GEOMETRY);

          if (cell->root_object == NULL) {
            cell->root_object = new BarrList();
            cell->object_type = BARR_LIST;
          }

          ((BarrList*)(cell->root_object))->list.push_back(new_obj);
        }
      }
    }
  }
  bool intersect(const Ray& r, float t_min, float t_max, Intersection& isect) {
    glm::vec3 m = bbox.minimum;

    glm::vec3 cell_step = (bbox.maximum - bbox.minimum);
    cell_step.x /= x_divs;
    cell_step.y /= y_divs;
    cell_step.z /= z_divs;

    float t0;
    bool hit = bbox->hit_t(r, t_min, t_max, t0);
    if (!hit)
      return false;

    glm::vec3 pos = r.at(t0);



    return false;
  }

};
*/
