#ifndef CGL_STATICSCENE_MEDIA_H
#define CGL_STATICSCENE_MEDIA_H

#include "object.h"
#include "primitive.h"

namespace CGL { namespace StaticScene {

/**
 * Participating media primitive.
 * Want to represent grid of densities.
 */
 class Media : public Primitive {
 public:
   /**
    * Constructor.
    */
  Media(const Grid* grid, size_t x, size_t y, size_t z);

  /**
   * Get the world space bounding box of the participating media.
   * \return world space bounding box of the participating media.
   */
 BBox get_bbox() const;

 bool is_medium() const {
   return true;
 }

 float medium_dist(const Ray& r) const;

 float intersection(const Ray& r) const;


 }

}}
