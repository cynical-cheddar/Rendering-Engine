#include "Wireframe.h"


#define WIDTH 320
#define HEIGHT 240

pixelMatrixDictionaryEntry::pixelMatrixDictionaryEntry(std::vector<std::vector<uint32_t>>  pixelMatrix, std::string name):
	pixelMatrix(std::move(pixelMatrix)),
    name(std::move(name)){}


glm::vec3 transformForward = glm::vec3(0,0,1);
World world;
std::vector<Material> materials;
std::vector<PointLight> pointLights;
TextureMap environmentMap = TextureMap("week5Files/models/spacebox3.ppm");
std::vector<std::vector<uint32_t>>  environmentPixelMatrix;

glm::vec3 camPos;
glm::vec3 origin;
bool enableOrbit = false;
glm::mat3 rotationDefault = glm::mat3(
		1.0, 0.0, 0.0,
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0
	);

std::vector<pixelMatrixDictionaryEntry> pixelMatrixDictionary;






Camera cam = Camera(camPos, rotationDefault, 4.0);

enum moveKey{
	w,
	a,
	s,
	d,
	q, //rotate anti clockwise about model pos
	e, // rotate clockwise around model pos
	up,
	down,
	left,
	right,
	o,
	p,
	l,
	i,
	j,
	r,
	f
};

enum RenderMode{
	wireframe,
	rasterise,
	raytrace
};
RenderMode rendermode;



std::vector<int> getIndexesFromFaceVector(std::vector<std::string> line){
	std::vector<std::string> faceVectors1 = split(line[1], '/');
	std::vector<std::string> faceVectors2 = split(line[2], '/');
	std::vector<std::string> faceVectors3 = split(line[3], '/');
	std::vector<int> indexes;
	indexes.push_back(stoi(faceVectors1.at(0)));
	indexes.push_back(stoi(faceVectors2.at(0)));
	indexes.push_back(stoi(faceVectors3.at(0)));
	return indexes;
}

struct StringColourPair{
	std::string colourString;
	std::string textureString;
	bool textured;
	Colour colour;
};

struct VerticesTexturePointPair{
	glm::vec3 faceVector;
	TexturePoint texVector;
};

struct TrianglesAndPalette{
	std::vector<ModelTriangle> triangles;
	std::vector<StringColourPair> colourPalette;
};

Colour lookupColourFromPalette(std::vector<StringColourPair> &colourPalette, std::string &keyString){
	// loop through vector by keystring to find colour
	for(StringColourPair colourPair:colourPalette){
		if(colourPair.colourString == keyString){
			return colourPair.colour;
		}
	}
	Colour nullCol;
	return nullCol;
}

StringColourPair lookupStringColourPairFromPalette(std::vector<StringColourPair> &colourPalette, std::string &keyString){
	for(StringColourPair colourPair:colourPalette){
		if(colourPair.colourString == keyString){
			return colourPair;
		}
	}
	
	return colourPalette[0];
}


/*
TrianglesAndPalette loadOBJFileLegacy(float scaleFactor){
	TrianglesAndPalette returnStruct;
	std::vector<ModelTriangle> triangles;
	std::vector<VerticesTexturePointPair> vertices;
	std::vector<StringColourPair> colourPalette;
	// load file
	//src/cornell-box.obj
	//week5Files/models/textured-cornell-box.obj
	std::ifstream objFile("week5Files/models/textured-cornell-box.obj");
	std::ifstream objFileSecondPass("week5Files/models/textured-cornell-box.obj");
	std::ifstream mtlFile("week5Files/models/textured-cornell-box.mtl");
	// get a list of all vectors in the OBJ, starting from 1
	bool writingObject = false;
	// first, read the .mtl to define a colour pallate

	std::string currentColourName;
	while(!mtlFile.eof()){
		std::string line;
		getline(mtlFile,line);
		std::vector<std::string> materialLine = split(line, ' ');
		
		if (materialLine.at(0) == "newmtl"){
			currentColourName = materialLine.at(1);
		}
		if (materialLine.at(0) == "Kd"){
			// read a colour and add it to our colour Paletta
			Colour newCol;
			newCol.red = stof(materialLine.at(1))*255; 
			newCol.green = stof(materialLine.at(2))*255; 
			newCol.blue = stof(materialLine.at(3))*255; 
			StringColourPair newStringColourPair;
			newStringColourPair.colourString = currentColourName;
			newCol.name = currentColourName;
			newStringColourPair.colour = newCol;
			colourPalette.push_back(newStringColourPair);
			//std::cout << "Face Colour " <<newStringColourPair.colour <<std::endl;
		}
		if (materialLine.at(0) == "map_Kd"){
			std::cout << "map_Kd " << materialLine.at(1) <<std::endl;
			std::cout << "map_Kd " << currentColourName <<std::endl;
			colourPalette.back().textureString = materialLine.at(1);
			colourPalette.back().textured = true;
		}
	}
	while(!objFile.eof()){
		std::string line;
		getline(objFile,line);
		// firstly, we need to get a vector of all of the vertices
		
		if(line[0]=='v'){
			// if the line describes a vector, read it in and define it
			std::vector<std::string> vectors = split(line, ' ');
			// for each vector, turn it into a float
			float x = stof(vectors[1]);
			float y = stof(vectors[2]);
			float z = stof(vectors[3]);
			VerticesTexturePointPair newVector;
			//glm::vec3 newVector;
			newVector.faceVector.x = x*scaleFactor; newVector.faceVector.y=y*scaleFactor; newVector.faceVector.z = z*scaleFactor;
			vertices.push_back(newVector);
		}
		if(line[0] == 'v' && line[1] == 't'){
			std::vector<std::string> vectors = split(line, ' ');
			// for each vector, turn it into a float
			float x = stof(vectors[1]);
			float y = stof(vectors[2]);

			vertices.back().texVector.x = x*255;
			vertices.back().texVector.y = y*255;

		}
		
	}
	Colour currentColour;
	
	while(!objFileSecondPass.eof()){
		std::string line;
		getline(objFileSecondPass,line);
		// split line, if it contains usemtl, then define a material to use for the following triangles
		std::vector<std::string> materialLine = split(line, ' ');
		
		// get the current colour from our pallete by string
		if (materialLine.at(0) == "usemtl"){
			std::string currentColourString = materialLine.at(1);
			//std::cout << currentColourString<<std::endl;
			currentColour = lookupColourFromPalette(colourPalette, currentColourString);
			//std::cout << currentColour<<std::endl;
		}
		if(line[0]=='f'){
			// if we are describing a face, then define a new triangle
			ModelTriangle newTriangle;
			
			// get the vertex indexes from the line
			std::vector<std::string> vectorIndexes = split(line, ' ');
			std::vector<int> indexes = getIndexesFromFaceVector(vectorIndexes); 
			// look up the verticies from our list
			glm::vec3 p1 = vertices.at(indexes.at(0)-1).faceVector;
			glm::vec3 p2 = vertices.at(indexes.at(1)-1).faceVector;
			glm::vec3 p3 = vertices.at(indexes.at(2)-1).faceVector;
			TexturePoint t1 = vertices.at(indexes.at(0)-1).texVector;
			TexturePoint t2 = vertices.at(indexes.at(1)-1).texVector;
			TexturePoint t3 = vertices.at(indexes.at(2)-1).texVector;
			
			std::array<glm::vec3,3> vecs;
			std::array<TexturePoint,3> vecsT;
			vecs[0] = p1;
			vecs[1] = p2;
			vecs[2] = p3;
			vecsT[0] = t1;
			vecsT[1] = t2;
			vecsT[2] = t3;
			newTriangle.vertices = vecs;
			newTriangle.texturePoints = vecsT;




			// get the current colour to paint the 
			newTriangle.colour = currentColour;
			triangles.push_back(newTriangle);
			//std::cout << newTriangle << std::endl;
			
		}
	}
	//std::cout << triangles.at(0) << std::endl;
	//std::cout << triangles.at(1) << std::endl;
	//std::cout << triangles.at(2) << std::endl;
	//	std::cout << triangles.at(3) << std::endl;
	returnStruct.triangles = triangles;
	returnStruct.colourPalette = colourPalette;
	return returnStruct;

}
*/

glm::vec3 localPointToScreenPos(glm::vec3 &point, float focalLength, float scale){
	float u = (focalLength * (point.x*-scale/point.z)) + (WIDTH/2);
	float v = (focalLength * (point.y*scale/point.z)) + (HEIGHT/2);





	glm::vec3 screenPos;
	screenPos.x = u;
	screenPos.y = v;
	screenPos.z = point.z;
	return screenPos;
}

std::vector<glm::vec3> modelTriangleToScreenCoordinates(ModelTriangle &triangle, glm::vec3 &camPos, float focalLength){
	// convert the world space coordinates to camera local space
	float magnificationFactor = 100;
	glm::vec3 p1_local = triangle.vertices.at(0) - camPos;
	glm::vec3 p2_local = triangle.vertices.at(1) - camPos;
	glm::vec3 p3_local = triangle.vertices.at(2) - camPos;

	// now that we have the local space coordianates of the triangle, we should calculate the 2d space
	glm::vec3 p1_screenPos = localPointToScreenPos(p1_local, focalLength, magnificationFactor);
	glm::vec3 p2_screenPos = localPointToScreenPos(p2_local, focalLength, magnificationFactor);
	glm::vec3 p3_screenPos = localPointToScreenPos(p3_local, focalLength, magnificationFactor);

	// package points up as a vector
	std::vector<glm::vec3> coordinateSelection;
	coordinateSelection.push_back(p1_screenPos);
	coordinateSelection.push_back(p2_screenPos);
	coordinateSelection.push_back(p3_screenPos);
	return coordinateSelection;
}
CanvasPoint vectorToCanvasPoint(glm::vec3 point){
	CanvasPoint p;
	p.x = point.x;
	p.y = point.y;
	p.depth = point.z;
	return p;
}

ModelTriangle adjustTriangleRotationPerspective(Camera camera, ModelTriangle triangle){
	ModelTriangle newModelTriangle;
	
	newModelTriangle = triangle;
	int i = 0;
	for( glm::vec3 vertex : newModelTriangle.vertices){
		glm::vec3 vertexToCam = vertex - camera.position;
		glm::vec3 adjustedVector  =  camera.rotation * vertexToCam;
		newModelTriangle.vertices[i] = adjustedVector + camera.position;
		i++;
	}
	return newModelTriangle;
	
}

Camera LookAt(Camera camLocal, glm::vec3 targetPos){
	Camera newCam = camLocal;
	glm::vec3 newForward =  camLocal.position - targetPos;

	newForward = glm::normalize(newForward);
	// new right is (0,1,0) x forward
	glm::vec3 newRight = glm::cross(newForward, glm::normalize(glm::vec3(0.0,1.0,0.0)));
	// new up is forward x right
	glm::vec3 newUp = glm::cross(newForward, newRight);
	
	glm::mat3 newRot = glm::mat3(
		newRight, -newUp, newForward
	);

	

	
	newCam.rotation = newRot;
	return newCam;
}


void projectObjToImagePlane(World world, std::vector<Material> materials , DrawingWindow &window, Camera &camera, RenderMode rendermode){
	
	std::vector<std::vector<float>> depthBuffer(WIDTH, std::vector<float>(HEIGHT, 0));
	
	glm::vec3 camPos = camera.position;

	float focalLength = camera.focalLength;
	
	

	int i = 0;
	// loop through each model triangle and place a pixel at the projected position
	for(WorldObject &worldObject : world.worldObjects){
		for(Object &object:worldObject.childObjects){
			for(ModelTriangle &tri : object.tris){
				tri.worldObjectId = i;
				ModelTriangle adjustedTri = adjustTriangleRotationPerspective(camera, tri);
				
				

				std::vector<glm::vec3> screenPosCoordinates = modelTriangleToScreenCoordinates(adjustedTri, camPos, focalLength);
				CanvasTriangle strokedTriangle;
				strokedTriangle.vertices[0] = vectorToCanvasPoint(screenPosCoordinates.at(0));
				strokedTriangle.vertices[1] = vectorToCanvasPoint(screenPosCoordinates.at(1));
				strokedTriangle.vertices[2] = vectorToCanvasPoint(screenPosCoordinates.at(2));
				
			//	Material mat = findMaterial(world.worldMaterials.at(i), object.material);
				Material mat = findMaterial(materials, object.material);
				
				if(mat.hasTexture){
				//	depthBuffer = drawRasterisedTriangleDepthBuffered(strokedTriangle, mat.colour, window, depthBuffer);
					// map the texture points

					strokedTriangle.vertices[0].texturePoint = tri.texturePoints[0];
					strokedTriangle.vertices[1].texturePoint = tri.texturePoints[1];
					strokedTriangle.vertices[2].texturePoint = tri.texturePoints[2];
				//	std::cout << "matname" << std::endl;
				//	std::cout << mat.name << std::endl;
					depthBuffer = drawDepthBufferedStrokedTriangle(strokedTriangle, mat.colour, window, depthBuffer);
					if(rendermode == RenderMode(rasterise)) depthBuffer = drawTexturedTriangleDepthBuffered(strokedTriangle, window, mat.texture, depthBuffer);
					//depthBuffer = drawRasterisedTriangleDepthBuffered(strokedTriangle, mat.colour, window, depthBuffer);
				}
				else if(!mat.hasTexture){
					depthBuffer = drawDepthBufferedStrokedTriangle(strokedTriangle, mat.colour, window, depthBuffer);
				//	drawRasterisedTriangle(strokedTriangle, mat.colour, window);
					if(rendermode == RenderMode(rasterise)) depthBuffer = drawRasterisedTriangleDepthBuffered(strokedTriangle, mat.colour, window, depthBuffer);
				}
			}
		}
		i++;
	}
}
		
glm::vec3 rotateCamAboutPositionX(glm::vec3 cameraPos, glm::vec3 rotationOrigin, float amt){
	glm::vec3 newPos = cameraPos;
	glm::mat3 rotationMatrix = glm::mat3(
		1.0, 0.0, 0.0,
		0.0, cos(amt), -sin(amt),
		0.0, sin(amt), cos(amt)
	);

	newPos = newPos * rotationMatrix;
	std::cout << "pos before: " << cameraPos.x <<"," << cameraPos.y << ","<< cameraPos.z << std::endl;
	std::cout << "pos after: " << newPos.x <<"," << newPos.y << ","<<newPos.z << std::endl;
	return newPos;
	
	
}
glm::vec3 rotateCamAboutPositionY(glm::vec3 cameraPos, glm::vec3 rotationOrigin, float amt){
	glm::vec3 newPos = cameraPos;
	glm::mat3 rotationMatrix = glm::mat3(
		cos(amt), 0.0, sin(amt),
		0.0, 1.0, 0,
		-sin(amt), 0.0, cos(amt)
	);
	
	newPos = (newPos-rotationOrigin) * rotationMatrix;
	std::cout << "pos before: " << cameraPos.x <<"," << cameraPos.y << ","<< cameraPos.z << std::endl;
	std::cout << "pos after: " << newPos.x <<"," << newPos.y << ","<<newPos.z << std::endl;
	return newPos;
	
	
}

Camera updateCameraRotationX(Camera cam, glm::vec3 artifactPos, float angle){
	Camera newCam = cam;
	
	glm::mat3 rotationMatrix(glm::vec3(1.0, 0.0, 0.0),
							 glm::vec3(0.0, cos(angle), sin(angle)),
							 glm::vec3(0.0, -sin(angle), cos(angle)));

	newCam.rotation = rotationMatrix * newCam.rotation;
	return newCam;
}
Camera updateCameraRotationY(Camera cam, glm::vec3 artifactPos, float angle){
	Camera newCam = cam;
	


	glm::mat3 rotationMatrix(glm::vec3(cos(angle), 0.0, -sin(angle)),
							 glm::vec3(0.0, 1.0, 0.0),
							 glm::vec3(sin(angle), 0.0, cos(angle)));

	newCam.rotation = rotationMatrix * newCam.rotation;
	return newCam;
}

void assignWorldObjectMaterials(WorldObject &worldObjects, std::vector<Material> &materials){
	for(Object &object:worldObjects.childObjects){
		for(ModelTriangle &tri : object.tris){
			
		
			Material mat = findMaterial(materials, object.material);
			tri.colour = mat.colour;
			tri.reflectivity = mat.reflectionValue;
			tri.opticalIndex = mat.opticalIndex;
			tri.refractiveMaterial = mat.refractive;
			
			
			
			if(mat.hasTexture){
				
				// assign texture to dictionary
				pixelMatrixDictionary.push_back(pixelMatrixDictionaryEntry(getPixelMatrix(mat.texture), mat.name));
			
			}
			if(mat.hasNormalMap){
				// assign nmap to dictionary
				std::cout << "we have a normal map" << std::endl;
				pixelMatrixDictionary.push_back(pixelMatrixDictionaryEntry(getPixelMatrixNormal(mat.normalMap), mat.name.append("_n")));
			}
		}
	}
}
// gives the amount of light we should reflect
float frensel(glm::vec3 &incidentDir, glm::vec3 &surfaceNormal, float n1, float n2){
	// cosi clamp
	float kr = 0.0f;
	float cosi = glm::dot(incidentDir, surfaceNormal);
	if(cosi<-1) cosi = -1;
	if(cosi>1) cosi = 1;
	if(cosi>0){
		float temp = n1;
		n1 = n2;
		n2 = temp;
	}
	float sint = (n1/n2) * sqrtf(std::max(0.0f, 1 - cosi*cosi));
	if(sint >= 1){
		kr = 1;
	}
	else{
		float cost = sqrtf(std::max(0.0f, 1 - sint * sint));
		cosi = fabsf(cosi);
		float Rs = ((n2 * cosi) - (n1 * cost)) / ((n2 * cosi) + (n1 * cost)); 
        float Rp = ((n1 * cosi) - (n2 * cost)) / ((n1 * cosi) + (n2 * cost)); 
        kr = (Rs * Rs + Rp * Rp) / 2.0f; 
	}
	return kr;
}
Colour uintToColour(uint32_t input){
	
	uint32_t redMask   = 255 << 16;
	uint32_t greenMask = 255 << 8;
	uint32_t blueMask  = 255;


	uint32_t redComponent   = (input & redMask) >> 16;
	uint32_t greenComponent = (input & greenMask) >> 8;
	uint32_t blueComponent =  (input & blueMask);
	Colour colour = Colour(redComponent, greenComponent, blueComponent);

	return colour;
}
glm::vec3 uintToVector(uint32_t input){
	
	uint32_t xMask   = 255 << 16;
	uint32_t yMask = 255 << 8;
	uint32_t zMask  = 255;


	uint32_t xComponent   = (input & xMask) >> 16;
	uint32_t yComponent = (input & yMask) >> 8;
	uint32_t zComponent =  (input & zMask);
	glm::vec3 vec = glm::vec3 (xComponent, yComponent, zComponent);
	vec = glm::normalize(vec);
	return vec;
}
std::vector<std::vector<uint32_t>> retrievePixelMatrix(std::string &key){
	for(pixelMatrixDictionaryEntry &tex : pixelMatrixDictionary){
		if(tex.name == key){
			return tex.pixelMatrix;
		}
	}
}

glm::vec3 interpolateNormalmapNormalraytrace(ModelTriangle &tri, Material &mat, BarycentricIntersectionCoordinates &barycentrics){
	glm::vec3 normal;
	// get the pixel matrix from our dictionary, remember to use the suffix "_n" so as to return the normal map and not the texture map
	std::vector<std::vector<uint32_t>>  pixelMatrix = retrievePixelMatrix(mat.name.append("_n"));
	//std::vector<uint32_t>  pixelMatrix = mat.texture.pixels;
	// interpolate texturepoint based on u and v
	float u = barycentrics.u;
	float v = barycentrics.v;

	glm::vec2 t0 = glm::vec2(tri.texturePoints.at(0).x, tri.texturePoints.at(0).y);
    glm::vec2 t1 = glm::vec2(tri.texturePoints.at(1).x, tri.texturePoints.at(1).y);
    glm::vec2 t2 = glm::vec2(tri.texturePoints.at(2).x, tri.texturePoints.at(2).y);

    glm::vec2 result;
  //  u = 1 - u;
  //  v = 1 - v;
    float w = 1 - v - u;
    // use barycentric coordinates
    result = (u)*t1 + (v)*t2 + (w)*t0;
   // result = glm::normalize(result);
	//col = uintToColour (mat.texture.pixels.at(round(result.x + result.y*mat.texture.width)));
	normal = uintToVector(pixelMatrix[result.x][result.y]);

	return normal;
}

Colour interpolateTextureColourRaytrace(ModelTriangle &tri, Material &mat, BarycentricIntersectionCoordinates &barycentrics){
	Colour col;
	std::vector<std::vector<uint32_t>>  pixelMatrix = retrievePixelMatrix(mat.name);
	//std::vector<uint32_t>  pixelMatrix = mat.texture.pixels;
	// interpolate texturepoint based on u and v
	float u = barycentrics.u;
	float v = barycentrics.v;

	glm::vec2 t0 = glm::vec2(tri.texturePoints.at(0).x, tri.texturePoints.at(0).y);
    glm::vec2 t1 = glm::vec2(tri.texturePoints.at(1).x, tri.texturePoints.at(1).y);
    glm::vec2 t2 = glm::vec2(tri.texturePoints.at(2).x, tri.texturePoints.at(2).y);

	//t0 = glm::normalize(t0);
	//t1 = glm::normalize(t1);
	//t2 = glm::normalize(t2);

    glm::vec2 result;
  //  u = 1 - u;
  //  v = 1 - v;
    float w = 1 - v - u;
    // use barycentric coordinates
    result = (u)*t1 + (v)*t2 + (w)*t0;
   // result = glm::normalize(result);
	//col = uintToColour (mat.texture.pixels.at(round(result.x + result.y*mat.texture.width)));
	if(result.x > mat.texture.width-1) result.x = mat.texture.width -2;
	if(result.y > mat.texture.height-1) result.y = mat.texture.height -2;
	col = uintToColour(pixelMatrix[result.x][result.y]);

	return col;
}


Colour fireRayTrace(DrawingWindow &window, Camera camera, World &world, std::vector<Material> &materials, std::vector<PointLight> &pointLights, glm::vec3 &rayDir ,int x, int y, bool reflectedRay, float currentOpticalDensity, int curPass){
			glm::vec3 additiveNormal = glm::vec3(0,0,0);
			CanvasPoint intersectionCanvasPoint;

			intersectionCanvasPoint.x = x;
			intersectionCanvasPoint.y = HEIGHT-y -1;

			BarycentricIntersectionCoordinates barycentricIntersectionCoordinates;
			RayTriangleIntersection intersection = getClosestIntersection(camera, rayDir, world, barycentricIntersectionCoordinates);
			if(!intersection.isNull)
			{
		    Material triangleMat = findMaterial(materials, intersection.intersectedTriangle.colour.name);
			Colour col = intersection.intersectedTriangle.colour;
			if(triangleMat.hasTexture) col = interpolateTextureColourRaytrace(intersection.intersectedTriangle, triangleMat, barycentricIntersectionCoordinates);
		//	if(triangleMat.hasNormalMap)std::cout<< "we got a map" << std::endl;
			// REFRACTIVE RAY STUFF------------------------------------------------------
			if(triangleMat.refractive && curPass <4){
				

				glm::vec3 surfaceNormal = calculateNormal(intersection.intersectedTriangle);
				surfaceNormal = glm::normalize(surfaceNormal);
				
				// interpolate via phong shading if available
				if(intersection.intersectedTriangle.hasNormals) surfaceNormal = interpolateVertexNormal(intersection.intersectedTriangle, barycentricIntersectionCoordinates.u, barycentricIntersectionCoordinates.v);
				
				glm::vec3 incidentDir = intersection.intersectionPoint - camera.position;
				incidentDir = glm::normalize(incidentDir);
				// cosi clamp
				float cosi = glm::dot(incidentDir, surfaceNormal);
				if(cosi<-1) cosi = -1;
				if(cosi>1) cosi = 1;


				

				float c1 = glm::dot(surfaceNormal, incidentDir);
				float n1 = currentOpticalDensity;
				float n2 = triangleMat.opticalIndex;
				// colour pixel by frensel amount and by white
				// gives fraction of the light of the pixel to be reflected
				
				// now that we have kr, 



				// we have exited a medium, set back to 1
				
				bool internal = false;
				bool totalInternalReflection = false;
				

				if(currentOpticalDensity == n2){
					 surfaceNormal = -surfaceNormal; internal = true; 
					 if(intersection.intersectedTriangle.hasNormals) n2 = 1.0f;
				}
				float kr = frensel(incidentDir, surfaceNormal, n1, n2);
				

				float n = n1/n2;

				float k = 1 - n*n*(1-cosi*cosi);
				if(k<0)totalInternalReflection = true;
				if(kr<0) totalInternalReflection = true;

				float c2 = sqrt(1.0f - (pow(n, 2) * (1-pow(c1, 2)) ));
				

				// calculate direction of refraction ray
				glm::vec3 refractionDir = n * incidentDir + (n*c1 - c2) * surfaceNormal;
				refractionDir = glm::normalize(refractionDir);
				// calculate direction of reflection ray
 				glm::vec3 refelectionDir = incidentDir - 2.0f*surfaceNormal* (glm::dot(incidentDir, surfaceNormal));
				refelectionDir = glm::normalize(refelectionDir);
				
				

				// now fire the new ray with n2's optical density
				curPass++;
				Colour refractCol;
				Camera tempcam = Camera(intersection.intersectionPoint, rotationDefault, camera.focalLength);
				if(!totalInternalReflection)refractCol = fireRayTrace(window, tempcam, world, materials, pointLights, refractionDir, x, y, false, n2, curPass);
				
				// now fire the reflection ray
				
				Colour reflectCol;
				tempcam = Camera(intersection.intersectionPoint + refelectionDir*0.01f, rotationDefault, camera.focalLength);
				reflectCol = fireRayTrace(window, tempcam, world, materials, pointLights, refelectionDir, x, y, true, triangleMat.opticalIndex, curPass);
				

				Colour hitColour = Colour(reflectCol.red*kr +refractCol.red*(1-kr), (reflectCol.green*kr +refractCol.green*(1-kr)), (reflectCol.blue*kr +refractCol.blue*(1-kr)));
				
				if(kr > 0.98) hitColour= reflectCol;
				if(kr < 0) hitColour = reflectCol;
				return hitColour;
				
				//hitColor += reflectionColor * kr + refractionColor * (1 - kr); 
			}


			// REFLECTIVE RAY STUFF------------------------------------------------------
			else if(triangleMat.reflective && curPass < 4){
				//std::cout << "REFLECTION" << std::endl;
				// if we have a reflective triangle, find a new colour by firing a new ray

				// find out how reflective we are by using a percentage of our original colour
				float reflectivity = triangleMat.reflectionValue;
				Colour partialCol = triangleMat.colour;



				glm::vec3 surfaceNormal = calculateNormal(intersection.intersectedTriangle);
				surfaceNormal = glm::normalize(surfaceNormal);
				bool environmentLighting = false;
				if(intersection.intersectedTriangle.hasNormals){
					// interpolate via phong shading
					environmentLighting = true;
					surfaceNormal = interpolateVertexNormal(intersection.intersectedTriangle, barycentricIntersectionCoordinates.u, barycentricIntersectionCoordinates.v);
				}

				glm::vec3 incidentDir = intersection.intersectionPoint - camera.position;
				incidentDir = glm::normalize(incidentDir);

				glm::vec3 refelectionDir = incidentDir - 2.0f*surfaceNormal* (glm::dot(incidentDir, surfaceNormal));
				Camera tempcam = Camera(intersection.intersectionPoint, rotationDefault, camera.focalLength);
			//	intersection = getClosestIntersection(tempcam, refelectionDir, world, barycentricIntersectionCoordinates);
				curPass++;
				return fireRayTrace(window, tempcam, world, materials, pointLights, refelectionDir, x, y, environmentLighting, currentOpticalDensity, curPass);
				//Colour col = intersection.intersectedTriangle.colour;
			}


			// RAY STUFF WHEN WE HIT A NOT SPECIAL SURFACE ----------------------
			else{

			// to do normal mapping, check if the material that we hit is normal mapped. if it is, then add the normal to wherever we hit
			// then normalise
			
			if(triangleMat.hasNormalMap){
				additiveNormal = interpolateNormalmapNormalraytrace(intersection.intersectedTriangle,triangleMat, barycentricIntersectionCoordinates);
			}
			
			// now get the world space of the intersection point we just traced
			// create a temporary camera object at its position
			curPass++;
			Camera tempCam = Camera( intersection.intersectionPoint, rotationDefault, camera.focalLength );

			// for each light source, get the camera to look at it
			for(PointLight &light : pointLights){
				tempCam = LookAt(tempCam, light.position);
				float intensitySum = 0;
				if(!light.soft){
					LightRayIntersection  lightRayIntersection = castLightRay(tempCam, world, light, intersection.intersectedTriangle, camera, barycentricIntersectionCoordinates, additiveNormal);
					if(lightRayIntersection.canSee) intensitySum += lightRayIntersection.relativeIntensity;
				}
				else{
					int n = 18;
					// create n lights in a sphere around the light point
					std::vector<glm::vec3> pointLightPositions = generatePointsInSphere(light.position, 0.3f, n);
					for(int i = 0; i < n; i++){
						glm::vec3 lightPos = pointLightPositions.at(i);
						PointLight miniLight;
						miniLight.position = lightPos;
						miniLight.range = light.range;
						miniLight.intensity = light.intensity / (float) n;
						LightRayIntersection  lightRayIntersection = castLightRay(tempCam, world, miniLight, intersection.intersectedTriangle, camera, barycentricIntersectionCoordinates, additiveNormal);
						if(lightRayIntersection.canSee) intensitySum += lightRayIntersection.relativeIntensity;

					}
				}
				if(intensitySum >1)intensitySum = 1;
				


				

				// for multiple lights, change this to a linear sum of weightsf
				if(intensitySum>AMBIENTLIGHTING){
					
						col.red = col.red*intensitySum;
						col.green = col.green*intensitySum;
						col.blue = col.blue*intensitySum;
				}
				else{
					col.red = col.red*AMBIENTLIGHTING;
					col.green = col.green*AMBIENTLIGHTING;
					col.blue = col.blue*AMBIENTLIGHTING;
				}
				
				//if we cannot see the light source, set the colour to be darker

			}
			return col;
			//drawPixel(intersectionCanvasPoint, col, window);
			}
			}
			// otherwise, draw the environment
			else if(reflectedRay){
				

				float m = 2.0f * sqrt( pow(rayDir.x, 2) + pow(rayDir.y, 2) + pow(rayDir.z +1, 2));

				float u = ((rayDir.x/m) + 0.5f)*environmentMap.width;
				float v = ((rayDir.y/m) + 0.5f)*environmentMap.height;

				// get colour from environment map
				if(u>environmentMap.width-1)u = environmentMap.width-1;
				if(v>environmentMap.height-1)v = environmentMap.height-1;
				if(v<0)v=0;
				if(u<0)u=0;
				uint32_t col = environmentPixelMatrix[u][v];
			/*	uint32_t redMask   = 255 << 16;
				uint32_t greenMask = 255 << 8;
				uint32_t blueMask  = 255;


				uint32_t redComponent   = (col & redMask);
				uint32_t greenComponent = (col & greenMask);
				uint32_t blueComponent =  (col & blueMask);*/
				Colour colour = uintToColour(col);

				return colour;
			//	window.setPixelColour(floor(x), floor(y), col);
				
				//drawPixel(intersectionCanvasPoint, col, window);
			}
			else{
				return Colour(0, 0, 0);
			
			}
			

			
		}			




void raytraceScene(DrawingWindow &window, Camera camera, World &world, std::vector<Material> &materials){
	window.clearPixels();
	
	// firstly, assign materials to all of our world objects
	for(WorldObject &worldObjects : world.worldObjects){
		assignWorldObjectMaterials(worldObjects, materials);
	}
	

	std::cout<<"raytracing"<<std::endl;
	std::cout<< world.worldObjects[0].childObjects[0].tris.size()<<std::endl;
	// now loop through each pixel of the scene and raytrace the scene
	for(float y = 0; y <HEIGHT; y++){
		//std::cout << y/(float)HEIGHT << "%" <<std::endl;
		for(float x = 0; x<WIDTH; x++){
			
			// calculate ray dir
			// focal length = f
			// width = 1
			// height = 1
			

			//check fraction along x with respect to -1 and 1
			float xAmt =  2*(x/WIDTH)-1;
			float yAmt = 2*(y/HEIGHT)-1;
			xAmt *= (4.0/3.0);

			glm::vec3 rayDir = glm::vec3(xAmt,yAmt,-camera.focalLength);
			rayDir = glm::normalize(rayDir);
			rayDir = rayDir * camera.rotation;
			Colour col = fireRayTrace(window, camera, world, materials, pointLights,rayDir, x, y, false, 1.0f, 0);
			CanvasPoint intersectionCanvasPoint;

			intersectionCanvasPoint.x = x;
			intersectionCanvasPoint.y = HEIGHT-y -1;
			drawPixel(intersectionCanvasPoint, col, window);
			
		}

	}
}
float lightY = 1;
float lightX = 0;

// ------------- WORLD STUFF ---------------------
void refreshWorld(World &worldCur){
	worldCur.worldObjects.clear();
	materials.clear();
	pointLights.clear();
}
void addObjectToWorld(std::string objPath, glm::vec3 position, float scale){
	std::ifstream objFile(objPath);
	WorldObject worldObject = readObjects(objFile, scale, position, "worldObject");
	world.worldObjects.push_back(worldObject);
}

void loadMaterial(std::string mtlPath){
	std::ifstream mtlFile(mtlPath);
	std::vector<Material> materialsLocal = readMaterials(mtlFile);
	materials.insert(std::end(materials), std::begin(materialsLocal), std::end(materialsLocal));
}
void addLightToWorld(PointLight pointLight){
	pointLights.push_back(pointLight);
}

// ------------ END WORLD STUFF

void renderScene(DrawingWindow &window){
	switch(rendermode){
		case wireframe:	projectObjToImagePlane(world, materials, window, cam, rendermode); break;
		case rasterise:	projectObjToImagePlane(world, materials, window, cam, rendermode); break;
		case raytrace:	raytraceScene(window, cam, world, materials); break;
	}
	//window.renderFrame();
}

void updateCameraPos(DrawingWindow &window, moveKey keyPressed){
	

	
	

//	PointLight p = PointLight(glm::vec3(lightX,lightY,1), 1.0, 9.0);
	refreshWorld(world);
  //  addLightToWorld (PointLight(glm::vec3(3,-1,7), 1, 18.0,  false, 1.0f));
	//addObjectToWorld("week5Files/models/box3.obj", glm::vec3(0,0,0), 0.7f);
	
	//addObjectToWorld("enterprise2.obj", glm::vec3(0,0.2,0), 1.0f);
	//loadMaterial("week5Files/models/textured-cornell-box.mtl");
	//loadMaterial("week5Files/models/logoMats.mtl");
	//loadMaterial("untitled.mtl");

	// test 4)
	//loadMaterial("week5Files/models/textured-cornell-box.mtl");
	//addLightToWorld (PointLight(glm::vec3(0,0.9,0), 1, 13.0,  false, 1.0f));
	//addLightToWorld (PointLight(glm::vec3(0,0.1,0), 1, 10.0,  false, 1.0f));
	//addObjectToWorld("week5Files/models/box3.obj",glm::vec3(0,0,0), 0.5f);
	//addObjectToWorld("week5Files/models/sphere.obj", glm::vec3(0,-0.9,0.75), 0.6f);

	//glm::vec3 enterprisePos = glm::vec3(0,0,38);
	//loadMaterial("week5Files/models/textured-cornell-box.mtl");
	//addLightToWorld ( PointLight( glm::vec3(0,0.7,0), 1, 40.0,  false, 1.0f));

	//addObjectToWorld("week5Files/models/box3.obj", glm::vec3(0,0,0), 0.6f);


	pointLights.clear();
	cam.position.z = 20;
	glm::vec3 enterpriseStart = glm::vec3(0,0,0);
	glm::vec3 enterprisePos = enterpriseStart;
	loadMaterial("week5Files/models/textured-cornell-box.mtl");
	glm::vec3 lightPos = (enterprisePos + glm::vec3(0,0.50f,0));
	addLightToWorld (PointLight(lightPos, 1, 10.0,  false, 1.0f));
	std::cout << "enterprisePos: " << enterprisePos.x <<","<< enterprisePos.y <<","<< enterprisePos.z << " , lightPos: " << lightPos.x <<","<< lightPos.y <<","<< lightPos.z << std::endl;
	addObjectToWorld("week5Files/models/box1.obj",enterprisePos, 0.5f);



		//addObjectToWorld("enterprise2.obj",enterprisePos, 1.0f);
	//addObjectToWorld("week5Files/models/logo.obj", glm::vec3(-10,-5,40), 0.025f);
	//cam = LookAt(cam, enterprisePos);


	//  2)
	//loadMaterial("week5Files/models/textured-cornell-box.mtl");
	//addLightToWorld (PointLight(glm::vec3(0,0.9,0), 1, 10.0,  true, 1.0f));
	//addObjectToWorld("week5Files/models/box1.obj",glm::vec3(0,0,0), 0.5f);
			//addObjectToWorld("week5Files/models/sphere2.obj", glm::vec3(0.5,0,0), 0.4f);
	
	window.clearPixels();

	Camera newCam = cam;

	switch(keyPressed){
		case w: newCam.position.z -=1; break;
		case a: newCam.position.x -=1; break;
		case s: newCam.position.z +=1; break;
		case d: newCam.position.x +=1; break;
		case up: newCam.position = rotateCamAboutPositionX(newCam.position, origin, -0.05); break;
		case down: newCam.position = rotateCamAboutPositionX(newCam.position, origin, 0.05); break;
		case left: newCam.position = rotateCamAboutPositionY(newCam.position, origin, 3.14/16); break;
		case right: newCam.position = rotateCamAboutPositionY(newCam.position, origin,- 3.14/16); break;
		case q: newCam = updateCameraRotationY(newCam, newCam.position, -3.14/16); break;
		case e: newCam = updateCameraRotationY(newCam, newCam.position, 3.14/16); break;
		case r: newCam = updateCameraRotationX(newCam, newCam.position, -3.14/16); break;
		case f: newCam = updateCameraRotationX(newCam, newCam.position, 3.14/16); break;
		case o: enableOrbit = !enableOrbit; break;
		case p: lightY += 0.1; break;
		case l: lightY -= 0.1; break;
		case i: lightX += 0.1; break;
		case j: lightX -= 0.1; break;
	}

	


	cam = newCam;

	
	std::cout << "rotation before LookAt" << std::endl;
	std::cout << glm::to_string(cam.rotation) << std::endl;
	glm::vec3 rotPoint = glm::vec3(0.0, 0.0, 0.0);
	
	std::cout << "rotation after LookAt" << std::endl;
	std::cout << glm::to_string(cam.rotation) << std::endl;
//	
	if(enableOrbit)cam = LookAt(cam, rotPoint);
//	std::cout << "checking first vertex normal" << std::endl;
	//std::cout << worldObject.childObjects[0].vertexNormals[0].x << worldObject.childObjects[0].vertexNormals[0].y <<  worldObject.childObjects[0].vertexNormals[0].z <<std::endl;
	switch(rendermode){
		case wireframe:	projectObjToImagePlane(world, materials, window, cam, rendermode); break;
		case rasterise:	projectObjToImagePlane(world, materials, window, cam, rendermode); break;
		case raytrace:	raytraceScene(window, cam, world, materials); break;
	}
	
}







	

	


void animateScene(float currentFrame, float fps, float duration, DrawingWindow &window, std::string sceneId){

	
	camPos.x = 0; camPos.y = 0; camPos.z = 28;


	cam.position = camPos;
	cam.rotation = rotationDefault;


	

	int stage = 0;
	float stage0End= 1.0f;
	float stage0Duration = 1.0f;

	float stage1Duration = 0.5f;
	float stage1End = 1.5f;

	float stage2Duration = 0.25f;
	float stage2End = 1.5f;

	
	float maxFrames = fps * duration;
	float t = currentFrame / maxFrames;
	float elapsedTime = (currentFrame / fps);

	if(elapsedTime>stage0End) stage++;
	if(elapsedTime>stage1End) stage++;
	if(elapsedTime>stage2End) stage++;

	glm::vec3 enterpriseStart = glm::vec3(0,0,0);

	glm::vec3 enterpriseEnd = glm::vec3(0,0,30);
	glm::vec3 enterpriseFinal = glm::vec3(0,0,58);

	glm::vec3 cornellStart = glm::vec3(-1,0,0);
	glm::vec3 cornellEnd = glm::vec3(3,0,25);
	glm::vec3 cornellEnd2 = glm::vec3(-3,0,25);
	glm::vec3 cornellMidBelow = glm::vec3(0,0,20);
	glm::vec3 cornellEndBelow = glm::vec3(0,-4,28);
	glm::vec3 camStart = glm::vec3(0,0,28);
	glm::vec3 camEnd = glm::vec3(-1,0,28);
	float startSphereScale = 0.25f;
	float endSphereScale = 1.2f;
	 RenderMode(wireframe);
		// change for rasterise to remove rotation
		if(sceneId == "0"){
		//	updateCameraPos(window, a);
			cam.position.x -=1; 
		}
		if(sceneId == "eR" || sceneId == "1" || sceneId == "2"){
		//	updateCameraPos(window, a);
			cam.position.x -=1; 
			//updateCameraPos(window, e);
			cam = updateCameraRotationY(cam, cam.position, 3.14/16);
		}
		if(sceneId == "eL" || sceneId == "4" || sceneId == "5"){
		//	updateCameraPos(window, d);
			cam.position.x +=1; 
			//updateCameraPos(window, q);
			cam = updateCameraRotationY(cam, cam.position, -3.14/16);
		}
		if(sceneId == "6"){
			cam.position.x -=2; 
		//	updateCameraPos(window, a);
		//	updateCameraPos(window, a);
		}
		

	// enterprise zoom
	if(currentFrame <= maxFrames && sceneId == "eR"){
		refreshWorld(world);
		window.clearPixels();
		if(stage == 0){
			float t1 = (elapsedTime/stage0Duration);
			
			if(t1 > 0.5) rendermode = RenderMode(wireframe);
			else rendermode = RenderMode(wireframe);

			addLightToWorld (PointLight(glm::vec3(camPos.x,0,camPos.z), 1, 50.0,  false, 1.0f));
			glm::vec3 enterprisePos = vectorLerp(enterpriseStart, enterpriseEnd, t1);
			loadMaterial("week5Files/models/textured-cornell-box.mtl");
			
			addObjectToWorld("enterprise2.obj",enterprisePos, 1.0f);
		}
	}
	if(currentFrame <= maxFrames && sceneId == "eL"){
		refreshWorld(world);
		window.clearPixels();
		if(stage == 0){
			float t1 = (elapsedTime/stage0Duration);
			
			if(t1 > 0.5) rendermode = RenderMode(raytrace);
			else rendermode = RenderMode(raytrace);

			addLightToWorld (PointLight(glm::vec3(1,0,28), 1, 30.0,  false, 1.0f));
			glm::vec3 enterprisePos = vectorLerp(enterpriseStart, enterpriseEnd, t1);
			loadMaterial("week5Files/models/textured-cornell-box.mtl");
			
			addObjectToWorld("enterprise2.obj",enterprisePos, 1.0f);
		}
	}
	// rasterised box
	else if(currentFrame <= maxFrames && sceneId == "0"){
		refreshWorld(world);
		window.clearPixels();
		if(stage == 0){
			float t1 = (elapsedTime/stage0Duration);
			
			if(t1 > 0.5) rendermode = RenderMode(rasterise);
			else rendermode = RenderMode(rasterise);

			
			glm::vec3 enterprisePos = vectorLerp(enterpriseStart+glm::vec3(-3,0,0), cornellEnd, t1);
			loadMaterial("week5Files/models/textured-cornell-box.mtl");
			
			addLightToWorld (PointLight(enterprisePos + glm::vec3(0,0.9,0), 1, 10.0,  true, 1.0f));
			addObjectToWorld("week5Files/models/box0.obj",enterprisePos, 0.5f);
		}
	
	}

	// mirror and soft lighting
	if(currentFrame <= maxFrames && sceneId == "1"){
		refreshWorld(world);
		window.clearPixels();
		if(stage == 0){
			float t1 = (elapsedTime/stage0Duration);
			
			if(t1 > 0.5) rendermode = RenderMode(raytrace);
			else rendermode = RenderMode(raytrace);

			pointLights.clear();
			glm::vec3 enterprisePos = vectorLerp(enterpriseStart, cornellEnd, t1);
			loadMaterial("week5Files/models/textured-cornell-box.mtl");
			glm::vec3 lightPos = (enterprisePos + glm::vec3(0,0.30f,0));
			addLightToWorld (PointLight(lightPos, 1, 10.0,  false, 1.0f));
			std::cout << "enterprisePos: " << enterprisePos.x <<","<< enterprisePos.y <<","<< enterprisePos.z << " , lightPos: " << lightPos.x <<","<< lightPos.y <<","<< lightPos.z << std::endl;
			addObjectToWorld("week5Files/models/box1.obj",enterprisePos, 0.5f);
		}
	}
	// normal mapping and frensel
	if(currentFrame <= maxFrames && sceneId == "2"){
		refreshWorld(world);
		window.clearPixels();
		if(stage == 0){
			float t1 = (elapsedTime/stage0Duration);
			
			if(t1 > 0.5) rendermode = RenderMode(wireframe);
			else rendermode = RenderMode(wireframe);
			glm::vec3 enterprisePos = vectorLerp(enterpriseStart, cornellEnd2, t1);
			loadMaterial("week5Files/models/textured-cornell-box.mtl");
			addLightToWorld ( PointLight( enterprisePos+ glm::vec3(0,0.9,0), 1, 13.0,  false, 1.0f));
			addObjectToWorld("week5Files/models/box2.obj", enterprisePos + glm::vec3(0,0,0), 0.5f);
			
		}
	}


	// DO THE ENTERPRISE LEFT (make separetely by changing the start directions)

	// mirror sphere and environment
	if(currentFrame <= maxFrames && sceneId == "4"){
		refreshWorld(world);
		window.clearPixels();
		if(stage == 0){
			float t1 = (elapsedTime/stage0Duration);
			
			if(t1 > 0.5) rendermode = RenderMode(wireframe);
			else rendermode = RenderMode(wireframe);

			
			glm::vec3 enterprisePos = vectorLerp(enterpriseStart, cornellEnd2, t1);
			loadMaterial("week5Files/models/textured-cornell-box.mtl");
			addLightToWorld ( PointLight( enterprisePos+ glm::vec3(0,0.9,0), 1, 13.0,  false, 1.0f));
			addObjectToWorld("week5Files/models/box3.obj", enterprisePos + glm::vec3(0,0,0), 0.5f);
			addObjectToWorld("week5Files/models/sphere2.obj", enterprisePos +  glm::vec3(0,-0.9,0.75), 0.6f);
		}
	}
	// glass sphere
	if(currentFrame <= maxFrames && sceneId == "5"){
		refreshWorld(world);
		window.clearPixels();
		if(stage == 0){
			float t1 = (elapsedTime/stage0Duration);
			
			if(t1 > 0.5) rendermode = RenderMode(wireframe);
			else rendermode = RenderMode(wireframe);

			
			glm::vec3 enterprisePos = vectorLerp(enterpriseStart, cornellEnd2, t1);
			loadMaterial("week5Files/models/textured-cornell-box.mtl");
			addLightToWorld ( PointLight( enterprisePos+ glm::vec3(0,0.9,0), 1, 13.0,  false, 1.0f));
			addObjectToWorld("week5Files/models/box3.obj", enterprisePos + glm::vec3(0,0,0), 0.5f);
			addObjectToWorld("week5Files/models/sphere2.obj", enterprisePos +  glm::vec3(0,-0.9,0.75), 0.6f);
		}
		
	}
	// this one will be three seconds long
	if(currentFrame <= maxFrames && sceneId == "6"){
		refreshWorld(world);
		window.clearPixels();
			float t1 = (elapsedTime/duration);
			
			if(t1 > 0.5) rendermode = RenderMode(wireframe);
			else rendermode = RenderMode(wireframe);

			
			glm::vec3 enterprisePos = vectorLerp(enterpriseStart, enterpriseFinal, t1);
			loadMaterial("week5Files/models/textured-cornell-box.mtl");
			addLightToWorld ( PointLight( glm::vec3(camPos.x*2, camPos.y, camPos.z), 1, 90.0,  false, 1.0f));
			addObjectToWorld("enterprise2.obj",enterprisePos, 1.0f);
			addObjectToWorld("week5Files/models/logo.obj", glm::vec3(-10,-5,62), 0.025f);
			cam = LookAt(cam, enterprisePos);
		
	}


	//addObjectToWorld("week5Files/models/sphere2.obj", glm::vec3(0.5,0,0), 0.4f);

	// rasterised box
	
	if(currentFrame<maxFrames){
		renderScene(window);
		window.renderFrame();
		std::string filename = sceneId;
		filename.append("_");
		
		filename.append(std::to_string(std::round(currentFrame)));
		filename.append(".bmp");
		window.saveBMP(filename);
		
	}
}
	


void draw(DrawingWindow &window) {
	//window.clearPixels();
	
}

void orbit(DrawingWindow &window){
	
	glm::vec3 rotPoint = glm::vec3(0.0, 0.0, 0.0);
	
	updateCameraPos(window, moveKey(left));
	
}

void update(DrawingWindow &window) {
	// Function for performing animation (shifting artifacts or moving the camera)
	// call orbit function
//	orbit(window);
	if(enableOrbit)orbit(window);
	
}
void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;

		if (event.key.keysym.sym == SDLK_w) updateCameraPos(window, moveKey(w));
		if (event.key.keysym.sym == SDLK_a) updateCameraPos(window, moveKey(a));
		if (event.key.keysym.sym == SDLK_s) updateCameraPos(window, moveKey(s));
		if (event.key.keysym.sym == SDLK_d) updateCameraPos(window, moveKey(d));
		if (event.key.keysym.sym == SDLK_UP) updateCameraPos(window, moveKey(up));
		if (event.key.keysym.sym == SDLK_DOWN) updateCameraPos(window, moveKey(down));
		if (event.key.keysym.sym == SDLK_LEFT) updateCameraPos(window, moveKey(left));
		if (event.key.keysym.sym == SDLK_RIGHT) updateCameraPos(window, moveKey(right));
		if (event.key.keysym.sym == SDLK_q) updateCameraPos(window, moveKey(q));
		if (event.key.keysym.sym == SDLK_e) updateCameraPos(window, moveKey(e));
		if (event.key.keysym.sym == SDLK_o) updateCameraPos(window, moveKey(o));
		if (event.key.keysym.sym == SDLK_p) updateCameraPos(window, moveKey(p));
		if (event.key.keysym.sym == SDLK_i) updateCameraPos(window, moveKey(i));
		if (event.key.keysym.sym == SDLK_j) updateCameraPos(window, moveKey(j));
		if (event.key.keysym.sym == SDLK_l) updateCameraPos(window, moveKey(l));
		if (event.key.keysym.sym == SDLK_r) updateCameraPos(window, moveKey(r));
		if (event.key.keysym.sym == SDLK_f) updateCameraPos(window, moveKey(f));
		if (event.key.keysym.sym == SDLK_1) rendermode = RenderMode(wireframe);
		if (event.key.keysym.sym == SDLK_2) rendermode = RenderMode(rasterise);
		if (event.key.keysym.sym == SDLK_3) rendermode = RenderMode(raytrace);
	} else if (event.type == SDL_MOUSEBUTTONDOWN) window.savePPM("output.ppm");
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	
	origin.x = 0; origin.y = 0; origin.z = 0;
	camPos.x = 0; camPos.y = 0; camPos.z = 28;


	cam.position = camPos;
	cam.rotation = rotationDefault;

	
	
	//std::ifstream mtlFile("week5Files/models/textured-cornell-box.mtl");
	
	//std::vector<Material> materials = readMaterials(mtlFile);

	//std::ifstream objFile("week5Files/models/textured-cornell-box.obj");
	//std::ifstream objFile("week5Files/models/sphere.obj");
	
	//WorldObject objects = readObjects(objFile, 0.3, origin, "cornell box");
	
	environmentPixelMatrix = getPixelMatrix(environmentMap);

	rendermode = RenderMode(wireframe);
	//updateCameraPos(window, moveKey(s));
	
	
	
	bool animate = true;
	if(animate == false){
		renderScene(window);
	} 

	//projectObjToImagePlane(objects, materials, window, cam);
	
	float currentFrame = 0;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		update(window);
		draw(window);
		if(animate){
		//	animateScene(currentFrame, 30.0f, 3.0f, window, "6");
		//	animateScene(currentFrame, 30.0f, 1.0f, window, "5");
		//	animateScene(currentFrame, 30.0f, 1.0f, window, "4");
		//	animateScene(currentFrame, 30.0f, 1.0f, window, "3");
		//	animateScene(currentFrame, 30.0f, 1.0f, window, "2");
			animateScene(currentFrame, 3.0f, 1.0f, window, "1");
		//	animateScene(currentFrame, 30.0f, 1.0f, window, "0"); [x]
		//	animateScene(currentFrame, 30.0f, 1.0f, window, "eR");
			//animateScene(currentFrame, 																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																					5.0f, 1.0f, window, "eL");
		} 
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		//std::cout << "CURRENT FRAME: " << currentFrame <<std::endl;
		window.renderFrame();
		currentFrame ++;
	}
	
}
