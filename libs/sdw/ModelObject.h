#pragma once

#include <string>
#include <vector>
#include <glm/glm.hpp>
#include "ModelTriangle.h"
#include "Utils.h"
#include <iostream>
#include <fstream>
#include "Colour.h"
#include "TextureMap.h"
#include "NormalMap.h"
#include "PointLight.h"

#define WIDTH 320
#define HEIGHT 240




struct Object {
    std::string name;
    std::string material;
    std::string parentName;
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> vertexNormals;
    std::vector<TexturePoint> texturePoints;
    std::vector<ModelTriangle> tris;
    
    
    Object();
    Object(std::string name);
    Object(std::string name, std::string material);
    Object(std::string name, std::string material, std::string parentName);
};

struct Material {
    bool hasNormalMap;
    bool reflective;
    bool refractive;
    float reflectionValue;
    float opticalIndex;
    std::string name;
    Colour colour;
    bool hasTexture;
    TextureMap texture;
    NormalMap normalMap;
    
    Material();
    Material(std::string name, Colour colour);
    Material(std::string name, Colour colour, TextureMap texture);
    Material(std::string name, Colour colour, NormalMap normalMap);
    Material(std::string name, Colour colour, float reflectionValue);
  //  Material(std::string name, Colour colour, float opticalIndex);
    Material(std::string name, Colour colour, float reflectionValue ,float opticalIndex);

    Material(std::string name, Colour colour , TextureMap texture , float reflectionValue);
    
//    Material(std::string name, Colour colour, TextureMap texture, float reflectionValue ,float opticalIndex);
};

struct WorldObject {
    std::string name;
    std::vector<Object> childObjects;
    glm::vec3 worldPosition;
    glm::mat3 worldRotation;
    WorldObject(std::string name,  std::vector<Object> childObjects, glm::vec3 worldPosition);
    WorldObject(std::string name,  std::vector<Object> childObjects, glm::vec3 worldPosition, glm::mat3 worldRotation);
};

struct World {
    std::vector<WorldObject> worldObjects;
    //std::vector<PointLight> lights;
   // std::vector<std::vector<Material>> worldMaterials;
    World();
  //  World(std::vector<WorldObject> worldObjects);
};

WorldObject readObjects(std::ifstream &objFile, float scale, glm::vec3 origin, std::string worldName);

std::vector<Material> readMaterials(std::ifstream &matFile);

const Material findMaterial(const std::vector<Material> &mats, std::string name);

void applyMaterials(std::vector<Object> &objs, std::vector<Material> &mats);