/*----------------------------------------------------------
* COSC363  Ray Tracer
*
*  The cylinder class
*  This is a subclass of SceneObject, and hence implements the
*  methods intersect() and normal().
-------------------------------------------------------------*/

#include "Cylinder.h"
#include <math.h>

/**
* Cylinder's intersection method.  The input is a ray. 
*/
float Cylinder::intersect(glm::vec3 p0, glm::vec3 dir)
{
    float a = dir.x * dir.x + dir.z * dir.z;
    float b = 2 * (dir.x * (p0.x - center.x) + dir.z * (p0.z - center.z));
    float c = (p0.x - center.x) * (p0.x - center.x) + (p0.z - center.z) * (p0.z - center.z) - radius * radius;

    float delta = b * b - 4 * a * c;

    if (delta < 0.001)
        return -1.0;

    float t1 = (-b - sqrt(delta)) / (2 * a);
    float t2 = (-b + sqrt(delta)) / (2 * a);

    float y1 = p0.y + t1 * dir.y;
    float y2 = p0.y + t2 * dir.y;

    if (y1 > center.y && y1 < center.y + height && t1 > 0)
        return t1;
    else if (y2 > center.y && y2 < center.y + height && t2 > 0)
        return t2;
    else
        return -1;
}

glm::vec3 Cylinder::normal(glm::vec3 p)
{
    glm::vec3 n(p.x - center.x, 0, p.z - center.z);
    n = glm::normalize(n);
    return n;
}

