#include "ModelTriangle.h"
#include <utility>
#include <strings.h>

ModelTriangle::ModelTriangle() = default;

ModelTriangle::ModelTriangle(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, Colour trigColour) :
		vertices({{v0, v1, v2}}), texturePoints(), colour(std::move(trigColour)), normal(), refractiveMaterial(false) {}



/*
ModelTriangle::ModelTriangle(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, Colour trigColour, std::string parentName) :
		vertices({{v0, v1, v2}}), texturePoints(), colour(std::move(trigColour)), normal(),  parentName(std::move(parentName)){}*/

std::ostream &operator<<(std::ostream &os, const ModelTriangle &triangle) {
	os << "(" << triangle.vertices[0].x << ", " << triangle.vertices[0].y << ", " << triangle.vertices[0].z << ")\n";
	os << "(" << triangle.vertices[1].x << ", " << triangle.vertices[1].y << ", " << triangle.vertices[1].z << ")\n";
	os << "(" << triangle.vertices[2].x << ", " << triangle.vertices[2].y << ", " << triangle.vertices[2].z << ")\n";
	return os;
}
