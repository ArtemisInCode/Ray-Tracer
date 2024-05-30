/*----------------------------------------------------------
* COSC363  Ray Tracer
*
*  The cylinder class
*  This is a subclass of SceneObject, and hence implements the
*  methods intersect() and normal().
-------------------------------------------------------------*/

#ifndef H_CYLINDER
#define H_CYLINDER

#include <glm/glm.hpp>
#include "SceneObject.h"

class Cylinder : public SceneObject
{
private:
    glm::vec3 baseCenter;
    float radius;
    float height;

public:
    Cylinder() {};
    Cylinder(glm::vec3 base, float r, float h) : baseCenter(base), radius(r), height(h) {}

    float intersect(glm::vec3 p0, glm::vec3 dir);
    glm::vec3 normal(glm::vec3 p);
};

#endif //!H_CYLINDER