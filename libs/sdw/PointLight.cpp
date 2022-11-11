#include "PointLight.h"

PointLight::PointLight() = default;

PointLight::PointLight(glm::vec3 position, float intensity, float range) :
    position(std::move(position)),
    intensity(std::move(intensity)),
    range(std::move(range)){}

PointLight::PointLight(glm::vec3 position, float intensity, float range, bool soft, float softRadius) :
    position(std::move(position)),
    intensity(std::move(intensity)),
    range(std::move(range)),
    soft(std::move(soft)),
    softRadius(std::move(softRadius)){}
    