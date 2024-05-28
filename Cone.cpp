/*----------------------------------------------------------

* COSC363  Ray Tracer

*

*  The cone class

*  This is a subclass of SceneObject, and hence implements the

*  methods intersect() and normal().

-------------------------------------------------------------*/

#include "Cone.h"
#include <math.h>

/**
 * Cone's intersection method. The input is a ray.
 */
float Cone::intersect(glm::vec3 p0, glm::vec3 dir)
{
    glm::vec3 vdif = p0 - center;
    float tanTheta = (radius / height) * (radius / height);

    float a = dir.x * dir.x + dir.z * dir.z - tanTheta * dir.y * dir.y;
    float b = 2 * (vdif.x * dir.x + vdif.z * dir.z - tanTheta * vdif.y * dir.y);
    float c = vdif.x * vdif.x + vdif.z * vdif.z - tanTheta * vdif.y * vdif.y;

    float delta = b * b - 4 * a * c;

    if (delta < 0.001) return -1.0;  // No intersection

    float t1 = (-b - sqrt(delta)) / (2 * a);
    float t2 = (-b + sqrt(delta)) / (2 * a);

    // Determine the intersection points
    glm::vec3 hit1 = p0 + t1 * dir;
    glm::vec3 hit2 = p0 + t2 * dir;

    // Check if the intersection points are within the bounds of the cone height
    if (t1 > 0 && hit1.y >= center.y && hit1.y <= center.y + height) {
        return t1;
    }

    if (t2 > 0 && hit2.y >= center.y && hit2.y <= center.y + height) {
        return t2;
    }

    return -1.0;  // No valid intersection within the cone bounds
}

/**
 * Returns the unit normal vector at a given point.
 * Assumption: The input point p lies on the cone.
 */
glm::vec3 Cone::normal(glm::vec3 p)
{
    glm::vec3 n = p - center;
    n.y = -(radius / height) * sqrt(n.x * n.x + n.z * n.z);
    n = glm::normalize(n);
    return n;
}
