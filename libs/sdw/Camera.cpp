#include "Camera.h"

Camera::Camera(glm::vec3 position, glm::mat3 rotation, float focalLength) :
    position(std::move(position)),
    rotation(std::move(rotation)),
    focalLength(std::move(focalLength)){}
    