#include "bvh.h"

#include "CGL/CGL.h"
#include "static_scene/triangle.h"

#include <iostream>
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
    return node;
  }

  vector<Primitive *> left = vector<Primitive *>();
  vector<Primitive *> right = vector<Primitive *>();

  double max_c = max(centroid_box.extent.x, max(centroid_box.extent.y, centroid_box.extent.z));
  Vector3D split = centroid_box.centroid();

  for (Primitive *p: prims) {
    if (max_c == centroid_box.extent.x) {
      if (p->get_bbox().centroid().x < split.x) left.push_back(p);
      else right.push_back(p);
    } else if (max_c == centroid_box.extent.y) {
      if (p->get_bbox().centroid().y < split.y) left.push_back(p);
      else right.push_back(p);
    } else { // max_c == z axis
      if (p->get_bbox().centroid().z < split.z) left.push_back(p);
      else right.push_back(p);
    }
  }

  if (left.empty()) {
    std::move(right.begin() + right.size() / 2, right.end(), left.begin());
    right.erase(right.begin() + right.size() / 2, right.end());
  } else if (right.empty()) {
    std::move(left.begin() + left.size() / 2, left.end(), right.begin());
    left.erase(left.begin() + left.size() / 2, left.end());
  }

  node->l = construct_bvh(left, max_leaf_size);
  node->r = construct_bvh(right, max_leaf_size);

  return node;
}


bool BVHAccel::intersect(const Ray& ray, BVHNode *node) const {
  // Part 2, task 3: replace this.
  // Take note that this function has a short-circuit that the
  // Intersection version cannot, since it returns as soon as it finds
  // a hit, it doesn't actually have to find the closest hit.
  if (node == NULL) return false;

  if (!node->bb.intersect(ray, ray.min_t, ray.max_t)) return false;

  if (node->isLeaf()) {
    for (Primitive *p: *(node->prims)) {
      total_isects++;
      if (!p->is_medium()) {
        if (p->intersect(ray)) {
          return true;
        }
      }
    }
    return false;
  }

  bool hit_left = false, hit_right = false;
  if (node->l != NULL) hit_left = intersect(ray, node->l);
  if (node->r != NULL) hit_right = intersect(ray, node->r);
  if (hit_left || hit_right) return true;

  return false;

}

bool BVHAccel::intersect(const Ray& ray, Intersection* i, BVHNode *node) const {
  // Part 2, task 3: replace this
  if (node == NULL) return false;

  if (!node->bb.intersect(ray, ray.min_t, ray.max_t)) return false;

  bool hit = false;
  if (node->isLeaf()) {
    for (Primitive *p: *(node->prims)) {
      total_isects++;
      if (!p->is_medium()) {
        if (p->intersect(ray, i)) {
          hit = true;
        }
      }
    }
    return hit;
  }

  bool hit_left = false, hit_right = false;
  if (node->l != NULL) hit_left = intersect(ray, i, node->l);
  if (node->r != NULL) hit_right = intersect(ray, i, node->r);
  if (hit_left || hit_right) return true;

  return false;

}

bool BVHAccel::intersect_medium(const Ray& ray, BVHNode *node, double &d) const {
  if (node == NULL) return false;

  if (!node->bb.intersect(ray, ray.min_t, ray.max_t)) return false;

  bool hit = false;

  if (node->isLeaf()) {
    for (Primitive *p: *(node->prims)) {
      total_isects++;
      if (p->is_medium()) {
        double dis = p->medium_dist(ray);
        if (dis > 0) {
          d = dis;
          return true;
        }
      }
    }
    return hit;
  }

  bool hit_left = false, hit_right = false;
  if (node->l != NULL) hit_left = intersect_medium(ray, node->l, d);
  if (node->r != NULL) hit_right = intersect_medium(ray, node->r, d);
  if (hit_left || hit_right) return true;

  return false;
}

}  // namespace StaticScene
}  // namespace CGL
