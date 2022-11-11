#include <iostream>
#include <fstream>
#include <vector>
#include "SDL.h"

#include "CanvasTriangle.h"
#include "TextureMap.h"
#include "DrawingWindow.h"
#include "Utils.h"
#include "fstream"
#include "vector"
#include "CanvasPoint.h"
#include "Colour.h"
#include "RayTriangleIntersection.h"
#include "Camera.h"
#include "ModelObject.h"
#include "ModelTriangle.h"
#include "PointLight.h"

#define AMBIENTLIGHTING 0.15f

struct BarycentricIntersectionCoordinates{

    float u;
    float v;
    float w;
    BarycentricIntersectionCoordinates();
};

glm::vec3 calculateNormal(ModelTriangle &tri);

struct LightRayIntersection {
    bool hit;
    bool canSee;
    float relativeIntensity;
	glm::vec3 intersectionPoint;
	float distanceFromLight;
    ModelTriangle intersectedTriangle;
	LightRayIntersection();
	LightRayIntersection(bool hit, bool canSee, float relativeIntensity,  const glm::vec3 &intersectionPoint, float distanceFromLight, const ModelTriangle intersectedTriangle);

};

glm::vec3 vectorLerp(glm::vec3 A, glm::vec3 B, float t);

glm::vec3 interpolateVertexNormal(ModelTriangle &tri, float u, float v);

RayTriangleIntersection getClosestIntersection(Camera &camera, glm::vec3 dir, World &world, BarycentricIntersectionCoordinates &barycentricIntersectionCoordinates);

LightRayIntersection castLightRay(Camera &camera, World &world, PointLight &pointLight, ModelTriangle &litTriangle, Camera &mainCamera, BarycentricIntersectionCoordinates &barycentricIntersectionCoordinates, glm::vec3 &additiveNormal);