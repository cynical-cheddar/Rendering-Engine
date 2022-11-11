#pragma once

#include <CanvasTriangle.h>
#include <TextureMap.h>
#include <NormalMap.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <CanvasPoint.h>
#include <Colour.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <Utils.h>
#include <triangle.h>
//#include <ModelTriangle.h>
#include <ModelObject.h>
#include <glm/gtx/string_cast.hpp>
#include <Camera.h>
#include <RayTracing.h>
#include <RayTriangleIntersection.h>
#include <PointLight.h>
#include <cmath>
#include <random>

struct pixelMatrixDictionaryEntry{
    std::vector<std::vector<uint32_t>>  pixelMatrix;
    std::string name;
    pixelMatrixDictionaryEntry(std::vector<std::vector<uint32_t>>  pixelMatrix, std::string name);
};