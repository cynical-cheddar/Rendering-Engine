#include "RayTracing.h"
#include "triangle.h"



struct Ray {
    glm::vec3 direction;
    glm::vec3 startPoint;
  //  Ray(glm::vec3 direction, glm::vec3 startPoint);
};

BarycentricIntersectionCoordinates::BarycentricIntersectionCoordinates() = default;

LightRayIntersection::LightRayIntersection()  = default;

glm::vec3 calculateNormal(ModelTriangle &tri){
    glm::vec3 v0 = tri.vertices[0] ;
    glm::vec3 v1 = tri.vertices[1] ;
    glm::vec3 v2 = tri.vertices[2];

    glm::vec3 normal = glm::cross((v1 - v0),(v2 - v0));
    return normal;
}

glm::vec3 hitPoint(ModelTriangle &tri, float u, float v, WorldObject &worldObject){
    glm::vec3 localTriangleCoordinates;
    glm::vec3 modelTriangleCoordinates;
    glm::vec3 worldCoordinates;
    glm::vec3 e0 = tri.vertices[1] - tri.vertices[0]; // corresponds to u
    glm::vec3 e1 = tri.vertices[2] - tri.vertices[0]; // corresponds to v

    localTriangleCoordinates = tri.vertices[0] + u*e0 + v*e1;

    // now we need to transpose triangle coordianates to world coordiantes, relative to v0
    // transpose from trianle space to v0, then from model space to world space

    worldCoordinates = localTriangleCoordinates + worldObject.worldPosition;
  
    return worldCoordinates;
}

glm::vec3 vectorLerp(glm::vec3 A, glm::vec3 B, float t){
    return A*(1.f-t) + B*t;
}

glm::vec3 interpolateVertexNormal(ModelTriangle &tri, float u, float v){
    glm::vec3 vn0 = glm::normalize(tri.vertexNormals[0]);
    glm::vec3 vn1 = glm::normalize(tri.vertexNormals[1]);
    glm::vec3 vn2 = glm::normalize(tri.vertexNormals[2]);

    glm::vec3 result;
  //  u = 1 - u;
  //  v = 1 - v;
    float w = 1 - v - u;
    // use barycentric coordinates
    result = (u)*vn1 + (v)*vn2 + (w)*vn0;
    result = glm::normalize(result);


    
    // return
    return result;
}

LightRayIntersection castLightRay(Camera &camera, World &world, PointLight &pointLight, ModelTriangle &litTriangle, Camera &mainCamera, BarycentricIntersectionCoordinates &barycentricIntersectionCoordinates, glm::vec3 &additiveNormal){
    bool canSee = false;
    LightRayIntersection intersection;
    float range = pointLight.range;
    // firstly, convert our list of objects into a list of triangles
    std::vector<ModelTriangle> tris;
    int i = 0;
    for(WorldObject &worldObject : world.worldObjects){
        for(Object &obj : worldObject.childObjects){
            for(ModelTriangle &tri: obj.tris){
                tri.worldObjectId = i;
                tris.push_back(tri);
            }
        }
        i++;
    }

  //  std::cout<< tris[0].vertexNormals[0].x << "," << tris[0].vertexNormals[0].y << "," << tris[0].vertexNormals[0].z << std::endl;
    
   // Ray ray = Ray(camera.rotation * dir, camera.position);
    Ray ray;
    //ray.direction = dir*camera.rotation;
  //  ray.direction = glm::normalize(ray.direction);
    ray.direction =  pointLight.position - camera.position;
    ray.startPoint = camera.position;
    
    float t = 0;
    float u = 0;
    float v = 0;



    bool hit = false;
    // iterate through triangles to see if our matrix equation can be solved
    float bestT = INFINITY;
    ModelTriangle ourTri;
    for(ModelTriangle tri : tris){
        //std::cout << "iterate through tris:" <<  tri << std::endl;
        glm::vec3 e0 = tri.vertices[1] - tri.vertices[0];
        glm::vec3 e1 = tri.vertices[2] - tri.vertices[0];
        glm::vec3 SPVector = camera.position - tri.vertices[0];
        glm::mat3 DEMatrix(-ray.direction, e0, e1);
        glm::vec3 tuv = glm::inverse(DEMatrix) * SPVector;
        
        t = tuv.x;
        u = tuv.y;
        v = tuv.z;
    
        if((t>0)  && ((u+v) <= 1) && (u>= 0) && (v >= 0) && (u<=1) && (v<= 1) && !tri.refractiveMaterial){

            if(t<bestT && t > 0.000001){
                    intersection.intersectedTriangle = tri;
                        
                        
                        bestT = t;

                        intersection.distanceFromLight = t;
                        
                        ourTri = tri;
                        // now we should get the 3d coordinates
                        intersection.intersectionPoint = hitPoint(tri, u, v, world.worldObjects.at(tri.worldObjectId));
                        hit = true;
                }
            }
        }

    

    

    float pointLightDistance = glm::distance(pointLight.position, ray.startPoint);

    bool reflected = false;
    float overrideIntensity = AMBIENTLIGHTING;
    if(bestT < pointLightDistance && bestT < 0.999){
        // something is in the way of the light
        canSee = false;
        intersection.relativeIntensity = 0;
        // check if the material we intersected is a mirror. if it is, then reflect the ray to find a new light ray intersection
        // get new ray dir
        if(intersection.intersectedTriangle.reflectivity > 0.9){
            /*
            glm::vec3 incidentDir = ray.direction;
            incidentDir = glm::normalize(incidentDir);

            glm::vec3 surfaceNormal = calculateNormal(litTriangle);
            surfaceNormal = glm::normalize(surfaceNormal);

            // recalculate surface normal with phong shading
            
            if(litTriangle.hasNormals) surfaceNormal = interpolateVertexNormal(litTriangle, barycentricIntersectionCoordinates.u, barycentricIntersectionCoordinates.v);
            surfaceNormal = surfaceNormal + additiveNormal;
            
            surfaceNormal = glm::normalize(surfaceNormal);

            glm::vec3 refelectionDir = incidentDir - 2.0f*surfaceNormal* (glm::dot(incidentDir, surfaceNormal));

            // return reflection
            Camera newCam = Camera(intersection.intersectionPoint, glm::mat3(0,0,0,0,0,0,0,0,0) , camera.focalLength);

            // cast a shadow ray from the intersection point in the reflection dir.
            // the brightness of that object is our pixel brightness
         
            

            LightRayIntersection in = castLightRay(newCam, world, pointLight, litTriangle, mainCamera, barycentricIntersectionCoordinates, additiveNormal);
            reflected = true;
            if(in.relativeIntensity <AMBIENTLIGHTING) in.relativeIntensity = AMBIENTLIGHTING;
            else{ overrideIntensity = in.relativeIntensity; canSee = true;}*/
        }
    }
    else if(pointLightDistance>range){
        canSee = false;
        intersection.relativeIntensity = 0;
    }
    
    else{
        canSee = true;
        // if we can see, then calculate the intensity fraction of the ray
        // use range and intensity of light
        float rayDist = 0;
        rayDist = pointLightDistance;
        // calculate fraction of light's range the ray has travelled
        float distanceFraction = ((range-rayDist)/range);
        // obey inverse square law
        float intensityMultiplier = distanceFraction*distanceFraction;

        // now do angle based intensity stuff
        glm::vec3 surfaceNormal = calculateNormal(litTriangle);
        surfaceNormal = glm::normalize(surfaceNormal);

        // recalculate surface normal with gouraud shading
        
        if(litTriangle.hasNormals) surfaceNormal = interpolateVertexNormal(litTriangle, barycentricIntersectionCoordinates.u, barycentricIntersectionCoordinates.v);
        surfaceNormal = surfaceNormal + additiveNormal;
        
        surfaceNormal = glm::normalize(surfaceNormal);
        // ------------------- INCIDENCE LIGHTING -------------------------
        // dot product of surface normal and vector-to-light
        glm::vec3 vectorToLight = pointLight.position - ray.startPoint;
        vectorToLight = glm::normalize(vectorToLight);

        float dotProduct = glm::dot(surfaceNormal, vectorToLight);
        if(dotProduct<0) dotProduct = 0;

        // ------------------- SPECULAR LIGHTING -------------------------
        // vector of light to pixel
        glm::vec3 incidentDir = ray.startPoint - pointLight.position;
        incidentDir = glm::normalize(incidentDir);
        // vector of pixel to maincam
        glm::vec3 viewDir = mainCamera.position - ray.startPoint;
        viewDir = glm::normalize(viewDir);

        //specular resolution

        int n = 256;

        glm::vec3 refelectionDir = incidentDir - 2.0f*surfaceNormal* (glm::dot(incidentDir, surfaceNormal));
        refelectionDir = glm::normalize(refelectionDir);

        float specularIntensity = (float) pow(glm::dot(viewDir, refelectionDir), n) / 3;

        intersection.relativeIntensity = intensityMultiplier * dotProduct + (specularIntensity*distanceFraction);
        //std::cout << "pointlight intensity " << pointLight.intensity << std::endl;
        intersection.relativeIntensity *= pointLight.intensity;
        
        
    }

    if(intersection.relativeIntensity > 1 * pointLight.intensity) intersection.relativeIntensity = 1 * pointLight.intensity;
    //f(intersection.relativeIntensity < AMBIENTLIGHTING) intersection.relativeIntensity = AMBIENTLIGHTING;
    if(reflected) intersection.relativeIntensity = overrideIntensity;



    intersection.canSee = canSee;
    
    return intersection;
}

RayTriangleIntersection getClosestIntersection(Camera &camera, glm::vec3 dir, World &world, BarycentricIntersectionCoordinates &barycentricIntersectionCoordinates){
    RayTriangleIntersection intersection;
    intersection.isNull = true;
    
    float range = 10000000.0;
    // firstly, convert our list of objects into a list of triangles
    std::vector<ModelTriangle> tris;
    int i = 0;
    for(WorldObject &worldObject : world.worldObjects){
        for(Object &obj : worldObject.childObjects){
            for(ModelTriangle &tri: obj.tris){
                tri.worldObjectId = i;
                tris.push_back(tri);
            }
        }
        i++;
    }
    //std::cout << "here" << std::endl;
    
    

    // Define the ray
    
   // Ray ray = Ray(camera.rotation * dir, camera.position);
    Ray ray;
   // ray.direction = dir * camera.rotation;
    ray.direction = dir;
    ray.direction = glm::normalize(ray.direction);
    ray.startPoint = camera.position;
    float t = 0;
    float u = 0;
    float v = 0;


    
     

    // iterate through triangles to see if our matrix equation can be solved
    float bestT = INFINITY;
    for(ModelTriangle &tri : tris){
        //std::cout << "iterate through tris:" <<  tri << std::endl;
        glm::vec3 e0 = tri.vertices[1] - tri.vertices[0];
        glm::vec3 e1 = tri.vertices[2] - tri.vertices[0];
        glm::vec3 SPVector = camera.position - tri.vertices[0];
        glm::mat3 DEMatrix(-ray.direction, e0, e1);
        glm::vec3 tuv = glm::inverse(DEMatrix) * SPVector;
        
        t = tuv.x;
        u = tuv.y;
        v = tuv.z;
    
        if(t>0 && t<range && (u+v) <= 1 && u>= 0 && v >= 0 && u<=1 && v<= 1 && t > 0.0001){
         //   std::cout << "valid!" << std::endl;
            if(t<bestT){
                bestT = t;
                intersection.distanceFromCamera = t;
                intersection.intersectedTriangle = tri;
                
                
                // now we should get the 3d coordinates
                intersection.intersectionPoint = hitPoint(tri, u, v, world.worldObjects.at(tri.worldObjectId));
                barycentricIntersectionCoordinates.u = u;
                barycentricIntersectionCoordinates.v = v;
                barycentricIntersectionCoordinates.w = 1.0f - u - v;
                intersection.isNull = false;
             //   std::cout << "hit!" << std::endl;
            }
        }

    }

   // std::cout << intersection << std::endl;


    
    return intersection;
}