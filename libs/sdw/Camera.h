#pragma once
#include <glm/glm.hpp>
#include "Utils.h"

struct Camera {
    glm::vec3 position;
    glm::mat3 rotation;
    float focalLength;
    Camera(glm::vec3 position, glm::mat3 rotation, float focalLength);
};