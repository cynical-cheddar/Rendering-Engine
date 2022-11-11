#pragma once

#include <string>
#include <vector>
#include <glm/glm.hpp>

std::vector<std::string> split(const std::string &line, char delimiter);

std::vector<glm::vec3> generatePointsInSphere(glm::vec3 origin, float radius, int n);
