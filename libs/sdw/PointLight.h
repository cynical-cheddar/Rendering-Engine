#pragma once
#include <glm/glm.hpp>
#include "Utils.h"

struct PointLight {
    glm::vec3 position;
    float intensity;
    float range;
    bool soft;
    float softRadius;
    PointLight();
    PointLight(glm::vec3 position, float intensity, float range);
    PointLight(glm::vec3 position, float intensity, float range, bool soft, float softRadius);
};