#pragma once

#include <glm/glm.hpp>
#include <string>
#include <array>
#include "Colour.h"
#include "TexturePoint.h"


struct ModelTriangle {
	std::array<glm::vec3, 3> vertices{};
	std::array<TexturePoint, 3> texturePoints{};
	std::array<glm::vec3, 3> vertexNormals{};
	Colour colour{};
	glm::vec3 normal{};
	bool hasNormals;
	int worldObjectId;
	float reflectivity;
	float opticalIndex;
	bool refractiveMaterial;
//	std::string parentName;

	ModelTriangle();
	ModelTriangle(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, Colour trigColour);
	ModelTriangle(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, Colour trigColour, float reflectivity, float opticalIndex);
//	ModelTriangle(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, Colour trigColour, std::string parentName) ;
	friend std::ostream &operator<<(std::ostream &os, const ModelTriangle &triangle);
};
