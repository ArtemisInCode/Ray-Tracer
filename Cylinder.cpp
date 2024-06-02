/*----------------------------------------------------------
* COSC363  Ray Tracer
*
*  The cylinder class
*  This is a subclass of SceneObject, and hence implements the
*  methods intersect() and normal().
-------------------------------------------------------------*/

#include "Cylinder.h"
#include <cmath>

float Cylinder::intersect(glm::vec3 p0, glm::vec3 dir) {
    glm::vec3 d = dir;
    glm::vec3 vdif = p0 - baseCenter;
    d.y = 0; // Ignoring y-axis for cylinder body intersection
    vdif.y = 0;

    float a = glm::dot(d, d);
    float b = 2 * glm::dot(vdif, d);
    float c = glm::dot(vdif, vdif) - radius * radius;
    float delta = b * b - 4 * a * c;

    if (delta < 0) return -1.0;

    float t1 = (-b - sqrt(delta)) / (2 * a);
    float t2 = (-b + sqrt(delta)) / (2 * a);

    float y1 = p0.y + t1 * dir.y;
    float y2 = p0.y + t2 * dir.y;

    bool bodyHit = false;
    float tBody = -1.0;

    if (y1 > baseCenter.y && y1 < baseCenter.y + height) {
        tBody = t1;
        bodyHit = true;
    }
    if (y2 > baseCenter.y && y2 < baseCenter.y + height) {
        if (tBody < 0 || t2 < tBody) {
            tBody = t2;
            bodyHit = true;
        }
    }

    // Check for intersection with cap
    float tCap = (baseCenter.y + height - p0.y) / dir.y;
    glm::vec3 pCap = p0 + tCap * dir;
    if (glm::distance(glm::vec2(pCap.x, pCap.z), glm::vec2(baseCenter.x, baseCenter.z)) <= radius && tCap > 0) {
        if (tBody < 0 || tCap < tBody) {
            return tCap;
        }
    }

    return bodyHit ? tBody : -1.0;
}

glm::vec3 Cylinder::normal(glm::vec3 p) {
    // If point is on cap, return the normal pointing up
    if (p.y >= baseCenter.y + height - 0.001) {
        return glm::vec3(0, 1, 0);
    }

    glm::vec3 n = p - baseCenter;
    n.y = 0;
    return glm::normalize(n);
}
