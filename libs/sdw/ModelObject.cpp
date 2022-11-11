#include "ModelObject.h"

glm::mat3 rotationDef = glm::mat3(
		1.0, 0.0, 0.0,
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0
	);



Object::Object() = default;
Material::Material() = default;


Object::Object(std::string name) :
    name(std::move(name)) {}

Object::Object(std::string name, std::string material) :
    name(std::move(name)),
    material(std::move(material)){}


Object::Object(std::string name, std::string material, std::string parentName) :
    name(std::move(name)),
    
    material(std::move(material)),
    parentName(std::move(parentName)){}

Material::Material(std::string name, Colour colour) :
    name(std::move(name)),
    colour(std::move(colour)),
    hasNormalMap(false),
    reflective(false),
    refractive(false),
    reflectionValue(0.0f),
    opticalIndex(1.0f),
    hasTexture(false) {}

Material::Material(std::string name, Colour colour, TextureMap texture) :
    name(std::move(name)),
    colour(std::move(colour)),
    hasNormalMap(false),
    hasTexture(true),
    reflective(false),
    refractive(false),
    reflectionValue(0.0f),
    opticalIndex(1.0f),
    texture(std::move(texture)) {}

Material::Material(std::string name, Colour colour, float reflectionValue) :
    name(std::move(name)),
    colour(std::move(colour)),
    hasNormalMap(false),
    hasTexture(false),
    reflective(true),
    refractive(false),
    opticalIndex(1.0f),
    reflectionValue(std::move(reflectionValue)) {}

Material::Material(std::string name, Colour colour, NormalMap normalMap) :
    name(std::move(name)),
    colour(std::move(colour)),
    hasNormalMap(true),
    hasTexture(false),
    reflective(false),
    refractive(false),
    opticalIndex(1.0f),
    normalMap(std::move(normalMap)){}

// change later to allow for partial reflection multiplier
Material::Material(std::string name, Colour colour, float reflectionValue, float opticalIndex) :
    name(std::move(name)),
    hasNormalMap(false),
    colour(std::move(colour)),
    hasTexture(false),
    reflective(false),
    refractive(true),
    opticalIndex(std::move(opticalIndex)),
    reflectionValue(std::move(reflectionValue)) {}

Material::Material(std::string name, Colour colour , TextureMap texture , float reflectionValue) :
    name(std::move(name)),
    colour(std::move(colour)),
    hasNormalMap(false),
    hasTexture(true),
    reflective(true),
    refractive(false),
    opticalIndex(0.0f),
    texture(std::move(texture)),
    reflectionValue(std::move(reflectionValue)) {}

// change reflection stuff later
/*
Material::Material(std::string name, Colour colour , TextureMap texture , float reflectionValue ,float opticalIndex) :
    name(std::move(name)),
    colour(std::move(colour)),
    hasTexture(true),
    reflective(false),
    refractive(true),
    texture(std::move(texture)),
    opticalIndex(std::move(opticalIndex)),
    reflectionValue(std::move(reflectionValue)) {}*/







WorldObject::WorldObject(std::string name, std::vector<Object> childObjects, glm::vec3 worldPosition):
    name(std::move(name)),
    childObjects(std::move(childObjects)),
    worldPosition(std::move(worldPosition)){}

WorldObject::WorldObject(std::string name, std::vector<Object> childObjects, glm::vec3 worldPosition, glm::mat3 worldRotation):
    name(std::move(name)),
    childObjects(std::move(childObjects)),
    worldPosition(std::move(worldPosition)),
    worldRotation(std::move(worldRotation)){}

World::World() = default;


WorldObject readObjects(std::ifstream &objFile, float scale, glm::vec3 origin, std::string worldName) {
    std::string line;
    std::vector<Object> objs = std::vector<Object>();
    
    std::cout<< "origin: " << origin.x << " " << origin.y << " " << origin.z << std::endl;

    bool first = true;
    bool textured = false;
    int shift = 1;
    int nShift = 1;
    int n = 0;
    Object cObj = Object("renderObject", "defaultMaterial", worldName);
    //cObj.tris.reserve(100);
    int texCount = 0;
    int nCount = 0;
    while (!objFile.eof()) {
        std::getline(objFile, line);

        if (line.substr(0, 2) == "o ") {
            
            if (!first) objs.push_back(cObj);
            else first = false;
            textured = false;
            shift += cObj.vertices.size();
            nShift += n;
            std::cout << nShift << std::endl;
            Object cObj1;
            cObj = cObj1;
            cObj.name = line.substr(2);
            texCount = 0;
            nCount = 0;
            n = 0;
            
        //    cObj = Object(line.substr(2));
            
        }
        else if (line.substr(0, 7) == "usemtl ") {
            cObj.material = line.substr(7);
        }
        else if(line.substr(0, 3) == "vn "){
            std::vector<std::string> coords = split(line.substr(3), ' ');
            cObj.vertexNormals.push_back(
                glm::vec3(
                    (float)std::atof(coords[0].data()),
                    (float)std::atof(coords[1].data()),
                    (float)std::atof(coords[2].data())
                )
            );
            nCount++;
            n ++;
        }
        else if (line.substr(0, 3) == "vt ") {
            texCount ++;
            textured = true;
            std::vector<std::string> tCoords = split(line.substr(3), ' ');
            TexturePoint texturePoint;
            // change to texturepoint height?
            texturePoint.x = (float)std::atof(tCoords[0].data()) * WIDTH;
            texturePoint.y = (float)std::atof(tCoords[1].data()) * HEIGHT;
            cObj.texturePoints.push_back(texturePoint);
        }
        else if (line.substr(0, 2) == "v ") {
            std::vector<std::string> coords = split(line.substr(2), ' ');
            cObj.vertices.push_back(
                glm::vec3(
                    (float)std::atof(coords[0].data()) * scale,
                    (float)std::atof(coords[1].data()) * scale,
                    (float)std::atof(coords[2].data()) * scale
                ) + origin
            );
        }
        
        else if (line.substr(0, 2) == "f ") {
            // check if the line contains //. If it does, then it also has surface normals we need to deal with
           /* bool faceHasVertexNormals = false;
            if(line.find("//") != std::string::npos){
                faceHasVertexNormals = true;
            }*/
            
            ModelTriangle tri;

            

            std::vector<std::string> verts = split(line.substr(2), ' ');

            // split by space

            // foreach element in string
            std::string vertex0 = "";
            std::string vertex1 = "";
            std::string vertex2 = "";
            std::string tex0 = "";
            std::string tex1 = "";
            std::string tex2 = "";
            std::string norm0 = "";
            std::string norm1 = "";
            std::string norm2 = "";

            std::string vertex = "";
            std::string tex = "";
            std::string norm = "";

            int i = 0;
            
            for(std::string vertexInfo : verts){
                // read vertex info value by value
                if(i<3){
                    vertex = "";
                    tex = "";
                    norm = "";
                    int vIndex = 0;
                    for(char c : vertexInfo){
                        if(c != '/'){
                            
                            if(vIndex == 0){
                                if(i == 0)vertex0.push_back(c);
                                if(i == 1)vertex1.push_back(c);
                                if(i == 2)vertex2.push_back(c);
                            }
                            if(vIndex == 1){
                                if(i == 0)tex0.push_back(c);
                                if(i == 1)tex1.push_back(c);
                                if(i == 2)tex2.push_back(c);
                            }
                            if(vIndex == 2){
                                if(i == 0)norm0.push_back(c);
                                if(i == 1)norm1.push_back(c);
                                if(i == 2)norm2.push_back(c);
                            }
                        }
                        if(c == '/')vIndex ++;
                    }
                }
           //     if(i == 0) vertex0 = vertex; tex0 = tex; norm0= norm;
           //     if(i == 1) vertex1 = vertex; tex1 = tex; norm1= norm;
           //     if(i == 2) vertex2 = vertex; tex2 = tex; norm2= norm;
                i++;
            }
            
           // std::cout << "vertex0 index:" << (std::stoi(vertex0) - 1) << std::endl;
            tri.vertices[0] = cObj.vertices.at(std::stoi(vertex0) - shift);
            tri.vertices[1] = cObj.vertices.at(std::stoi(vertex1) - shift);
            tri.vertices[2] = cObj.vertices.at(std::stoi(vertex2) - shift);
           //  std::cout << "done vert" << std::endl;
         //    std::cout << "done vertex" << std::endl;
            
            if(tex0.length()>0 && texCount > 0){
                tri.texturePoints[0] = cObj.texturePoints.at(std::stoi(tex0) - 1);
                tri.texturePoints[1] = cObj.texturePoints.at(std::stoi(tex1) - 1);
                tri.texturePoints[2] = cObj.texturePoints.at(std::stoi(tex2) - 1);
              //  std::cout << "here:" << (std::stoi(vertex0) - 1) << std::endl;
            }
            else{
                TexturePoint p;
                p.x = -1;
                p.y = -1;
                tri.texturePoints[0] = p;
                tri.texturePoints[1] = p;
                tri.texturePoints[2] = p;
            }
           // std::cout << "done tex" << std::endl;
            if(norm0.length()>0 && nCount > 0){
               // std::cout << "read vertex normals" << std::endl;
                tri.vertexNormals[0] = cObj.vertexNormals.at(std::stoi(norm0) - nShift);
                tri.vertexNormals[1] = cObj.vertexNormals.at(std::stoi(norm1) - nShift);
                tri.vertexNormals[2] = cObj.vertexNormals.at(std::stoi(norm2) - nShift);
               // std::cout<< "vertexNormal0: " << tri.vertexNormals[0].x << " "<< tri.vertexNormals[0].y << " " <<tri.vertexNormals[0].z <<std::endl;
                tri.hasNormals = true;
            }
            else{
                glm::vec3 def = glm::vec3(0,0,0);
                tri.vertexNormals[0] = def;
                tri.vertexNormals[1] = def;
                tri.vertexNormals[2] = def;
                tri.hasNormals = false;
            }
          //  std::cout << "done norm" << std::endl;
          //   std::cout << "done norm" << std::endl;

          //  std::cout << "norm index" << (std::stoi(norm0) - 1) << std::endl;

            //std::cout<< "vertex0: " << tri.vertices[0].x << " "<< tri.vertices[0].y << " " <<tri.vertices[0].z <<std::endl;
        //    std::cout<< tri.texturePoints[0]<<std::endl;
      //      std::cout<< tri.vertexNormals[0].x << " " <<tri.vertexNormals[0].y << " " <<tri.vertexNormals[0].z <<std::endl;
            
            //std::vector<string> elements = split(Operations::trim(line.substr(2)), *"");
            /*
            std::vector<std::string> verts = split(line.substr(2), '/');
            ModelTriangle tri;
            tri.vertices[0] = cObj.vertices.at(std::atoi(verts[0].data()) - shift);
            tri.vertices[1] = cObj.vertices.at(std::atoi(verts[1].substr(1).data()) - shift);
            tri.vertices[2] = cObj.vertices.at(std::atoi(verts[2].substr(1).data()) - shift);
            tri.parentName = cObj.name;
            if(textured){
                // get the texture point for the current vertex
                int vertexNum = std::atoi(verts[0].data()) - shift;
                int vertexNum2 = std::atoi(verts[1].substr(1).data()) - shift;
                int vertexNum3 = std::atoi(verts[2].substr(1).data()) - shift;
                
                tri.texturePoints[0] = cObj.texturePoints.at(vertexNum);
                tri.texturePoints[1] = cObj.texturePoints.at(vertexNum2);
                tri.texturePoints[2] = cObj.texturePoints.at(vertexNum3);
                
               // std::cout<<"tri.texturePoints[0] - " << tri.texturePoints[0] << std::endl;
              //  std::cout<<"tri.texturePoints[1] - " << tri.texturePoints[1] << std::endl;
              //  std::cout<<"tri.texturePoints[2] - " << tri.texturePoints[2] << std::endl;
            }
            if(faceHasVertexNormals){
                int vertexNum = std::atoi(verts[0].data()) - shift;
                int vertexNum2 = std::atoi(verts[1].substr(1).data()) - shift;
                int vertexNum3 = std::atoi(verts[2].substr(1).data()) - shift;
                
                tri.vertexNormals[0] = cObj.vertexNormals.at(vertexNum);
                tri.vertexNormals[1] = cObj.vertexNormals.at(vertexNum2);
                tri.vertexNormals[2] = cObj.vertexNormals.at(vertexNum3);
            }
            */
            
            
            cObj.tris.push_back(tri);
           // std::cout<< "here" <<std::endl;
        }
    }
    cObj.parentName = worldName;
    objs.push_back(cObj);
    WorldObject worldObject("worldObject", objs, origin);
    worldObject.worldRotation = rotationDef;
    worldObject.name = worldName;
    return worldObject;
}

std::vector<Material> readMaterials(std::ifstream &matFile) {
    std::string line;
    std::vector<Material> mats = std::vector<Material>();
    bool first = true;
    std::string name;
    Colour colour;
    TextureMap tMap;
    NormalMap nMap;
    bool textureRead = false;
    bool reflectionRead = false;
    bool refractionRead = false;
    bool normalRead = false;
    float reflectivity = 0.0f;
    float opticalIndex = 1.0f;

    // add a defaultMaterial
    Colour defaultColour;
    defaultColour.name = "defaultColour";
    defaultColour.red = 255;
    defaultColour.green = 255;
    defaultColour.blue = 255;
    Material defaultMaterial = Material("defaultMaterial", defaultColour);
    mats.push_back(defaultMaterial);



    while (!matFile.eof()) {
        std::getline(matFile, line);

        if (line.substr(0, 7) == "newmtl ") {
            if (!first){
                if(normalRead){
                    std::cout << "normalRead" << std::endl;
                    mats.push_back(Material(colour.name, colour, nMap));
                }
                if(refractionRead){
                    std::cout << "refractionRead" << std::endl;
                    mats.push_back(Material(colour.name, colour, reflectivity, opticalIndex));
                }

               else if(reflectionRead && textureRead){
                   mats.push_back(Material(colour.name, colour, tMap, reflectivity));
               }

               else if(reflectionRead && !textureRead){
                   mats.push_back(Material(colour.name, colour, reflectivity));
               }
               
               else if(textureRead && !reflectionRead){
                    mats.push_back(Material(colour.name, colour, tMap));
               }
               else{
                   mats.push_back(Material(colour.name, colour));
               }
                reflectionRead = false;
                textureRead = false;
                refractionRead = false;
                normalRead = false;
                reflectivity = 0.0f;
                opticalIndex = 0.0f;
            } 
            else first = false;

            colour = Colour();
            colour.name = line.substr(7);
        }
        else if (line.substr(0, 3) == "Kd ") {
            std::vector<std::string> rgb = split(line.substr(3), ' ');
            colour.red = (int)(std::atof(rgb[0].data()) * 255);
            colour.green = (int)(std::atof(rgb[1].data()) * 255);
            colour.blue = (int)(std::atof(rgb[2].data()) * 255);
        }
        else if (line.substr(0, 4) == "Ref ") {
            std::vector<std::string> rgb = split(line.substr(4), ' ');
            reflectionRead = true;
            reflectivity = (float)(std::atof(rgb[0].data()));
        }
        else if (line.substr(0, 3) == "Ni ") {
            std::vector<std::string> rgb = split(line.substr(3), ' ');
            refractionRead = true;
            opticalIndex = (float)(std::atof(rgb[0].data()));
        }
        else if (line.substr(0, 7) == "map_Kd ") {
            std::vector<std::string> rgb = split(line.substr(7), ' ');
            std::string tMapName = rgb[0].data();
            std::cout << "tmapname " << tMapName << std::endl;
            tMap = TextureMap(tMapName);
            textureRead = true;
        }
        else if (line.substr(0, 7) == "map_Ns ") {
            std::vector<std::string> rgb = split(line.substr(7), ' ');
            std::string nMapName = rgb[0].data();
           
            nMap = NormalMap(nMapName);
            normalRead = true;
            std::cout << "nmapname " << nMapName << std::endl;
        }

    }
    if(textureRead)mats.push_back(Material(colour.name, colour, tMap));
    else mats.push_back(Material(colour.name, colour));
    

    return mats;
}

const Material findMaterial(const std::vector<Material> &mats, std::string name) {
    for (const Material &m : mats) {
        if (m.name == name) {
            return m;
        }
    }

    throw "fuck";
}

void applyMaterials(std::vector<Object> &objs, std::vector<Material> &mats) {
    for (Object &obj : objs) {
        const Colour &colour = findMaterial(mats, obj.material).colour;

        for (ModelTriangle &tri : obj.tris) {
            tri.colour = colour;
        }
    }
}