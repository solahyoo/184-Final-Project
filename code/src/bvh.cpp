#include "bvh.h"

#include "CGL/CGL.h"
#include "static_scene/triangle.h"

#include <iostream>
#include <limits>
#include <stack>

using namespace std;

namespace CGL { namespace StaticScene {

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  root = construct_bvh(_primitives, max_leaf_size);

}

BVHAccel::~BVHAccel() {
  if (root) delete root;
}

BBox BVHAccel::get_bbox() const {
  return root->bb;
}

void BVHAccel::draw(BVHNode *node, const Color& c) const {
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims))
      p->draw(c);
  } else {
    draw(node->l, c);
    draw(node->r, c);
  }
}

void BVHAccel::drawOutline(BVHNode *node, const Color& c) const {
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims))
      p->drawOutline(c);
  } else {
    drawOutline(node->l, c);
    drawOutline(node->r, c);
  }
}

/* Returns the index of the maximum value of vect */
int arg_max(Vector3D vect) {
  int index = 0;
  double max = std::numeric_limits<int>::min();
  for (int i = 0; i < 3; i++) {
    if (vect[i] > max) {
      max = vect[i];
      index = i;
    }
  }
  return index;
}

BVHNode *BVHAccel::construct_bvh(const std::vector<Primitive*>& prims, size_t max_leaf_size) {

  // Part 2, Task 1:
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.


  BBox centroid_box, bbox;

  for (Primitive *p : prims) {
    BBox bb = p->get_bbox();
    bbox.expand(bb);
    Vector3D c = bb.centroid();
    centroid_box.expand(c);
  }

  BVHNode *node = new BVHNode(bbox);
  if (prims.size() <= max_leaf_size) {
    node->prims = new vector<Primitive *>(prims);
  } else {
    int max_index = arg_max(centroid_box.extent);
    Vector3D midpoints = centroid_box.centroid();
    double split_point = midpoints[max_index];
    std::vector<Primitive *> *left = new vector<Primitive *>();
    std::vector<Primitive *> *right = new vector<Primitive *>();

    for (Primitive *p : prims) {
      if (split_point > p->get_bbox().centroid()[max_index]) {
        left->push_back(p);
      } else {
        right->push_back(p);
      }
    }
    if (left->size() == 0 || right->size() == 0) {
      node->prims = new vector<Primitive *>(prims);
    } else{
      node->l = construct_bvh(*left, max_leaf_size);
      node->r = construct_bvh(*right, max_leaf_size);
    }
  }
  return node;
}


bool BVHAccel::intersect(const Ray& ray, BVHNode *node) const {
  // Part 2, task 3: replace this.
  // Take note that this function has a short-circuit that the
  // Intersection version cannot, since it returns as soon as it finds
  // a hit, it doesn't actually have to find the closest hit.
  BBox bbox = node->bb;
  double t0, t1;
  if (!bbox.intersect(ray, t0, t1))
    return false;
  if (t1 < ray.min_t || t0 > ray.max_t)
    return false;
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims)) {
      total_isects++;
      if (p->intersect(ray))
        return true;
    }
    return false;
  }
  bool hit1 = intersect(ray, node->l);
  bool hit2 = intersect(ray, node->r);
  return hit1 || hit2;

}

bool BVHAccel::intersect(const Ray& ray, Intersection* i, BVHNode *node) const {
  // Part 2, task 3: replace this
  bool hit = false;
  BBox bbox = node->bb;
  double t0, t1;
  if (!bbox.intersect(ray, t0, t1))
    return false;
  if (t1 < ray.min_t || t0 > ray.max_t)
    return false;
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims)) {
      total_isects++;
      if (p->intersect(ray, i))
        hit = true;
    }
  } else{
    bool hit1 = intersect(ray, i, node->l);
    bool hit2 = intersect(ray, i, node->r);
    hit = hit1 || hit2;
  }
  return hit;
}

}  // namespace StaticScene
}  // namespace CGL
