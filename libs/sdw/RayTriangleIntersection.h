#pragma once

#include <glm/glm.hpp>
#include <iostream>
#include "ModelTriangle.h"

struct RayTriangleIntersection {
	glm::vec3 intersectionPoint;
	float distanceFromCamera;
	ModelTriangle intersectedTriangle;
	size_t triangleIndex;
	bool isNull;

	RayTriangleIntersection();
	RayTriangleIntersection(const glm::vec3 &point, float distance, const ModelTriangle &triangle, size_t index);
	friend std::ostream &operator<<(std::ostream &os, const RayTriangleIntersection &intersection);
	RayTriangleIntersection(const glm::vec3 &point, float distance, const ModelTriangle &triangle, size_t index, bool isNull);
	friend std::ostream &operator<<(std::ostream &os, const RayTriangleIntersection &intersection);
};
