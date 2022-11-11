#include "CanvasTriangle.h"
#include "TextureMap.h"
#include "DrawingWindow.h"
#include "Utils.h"
#include "fstream"
#include "vector"
#include "CanvasPoint.h"
#include "Colour.h"
#include "triangle.h"


#define WIDTH 320
#define HEIGHT 240

void drawPixel(CanvasPoint &point, Colour &col, DrawingWindow &window){
	float red = col.red;
	float green = col.green;
	float blue = col.blue;
	uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
	window.setPixelColour(floor(point.x), floor(point.y), colour);
}

void drawLine(CanvasPoint &point1, CanvasPoint &point2, Colour &col ,DrawingWindow &window){
	float xDiff = abs(point1.x - point2.x);
	float yDiff = abs(point1.y - point2.y);
	
	float xDiffSigned = point1.x - point2.x;
	float yDiffSigned = point1.y - point2.y;

	float stepCount= 0;
	if(xDiff  > yDiff) stepCount = xDiff;
	else stepCount = yDiff;
	float xStepSize = xDiffSigned/stepCount;
	float yStepSize = yDiffSigned/stepCount;

	for (float i = 0.0; i<stepCount; i++){
		float x = point2.x + (xStepSize*i);
		float y = point2.y + (yStepSize*i);

		float red = col.red;
		float green = col.green;
		float blue = col.blue;
		uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
		window.setPixelColour(floor(x), floor(y), colour);
	}
}
std::vector<float> interpolateDepthPoints(CanvasPoint point1, CanvasPoint point2, float n){
	std::vector<float> result;


	point1.depth = 1/point1.depth;
	point2.depth = 1/point2.depth;

	float zGap = (point1.depth-point2.depth)/(n-1);
	for (int i = 0; i < n; i++){
		float zNew = -(point2.depth + zGap*i);

		result.push_back(zNew);
	}
	return result;
}
std::vector<std::vector<float>> drawLineDepthBuffered(CanvasPoint &point1, CanvasPoint &point2, Colour &col ,DrawingWindow &window, std::vector<std::vector<float>> depthBuffer, float lineMultiplier){
	std::vector<std::vector<float>> newDepthBuffer = depthBuffer;
	
	float xDiff = abs(point1.x - point2.x);
	float yDiff = abs(point1.y - point2.y);
	
	float xDiffSigned = point1.x - point2.x;
	float yDiffSigned = point1.y - point2.y;
	

	float stepCount= 0;
	if(xDiff  > yDiff) stepCount = xDiff;
	else stepCount = yDiff;
	float xStepSize = xDiffSigned/stepCount;
	float yStepSize = yDiffSigned/stepCount;

	// now get a set of depth values size stepCount interpolated between p1 and p2
	std::vector<float> depthPoints = interpolateDepthPoints(point1, point2, stepCount);
	

	for (float i = 0.0; i<stepCount; i++){
		float x = point2.x + (xStepSize*i);
		float y = point2.y + (yStepSize*i);

		float red = col.red;
		float green = col.green; 
		float blue = col.blue;
		uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
		// check if the (1/z value of the current depth is greater than the stored position in the buffer)
		if(x > 0 && x < WIDTH && y > 0 && y < HEIGHT){
			if(depthPoints.at(floor(i)) > newDepthBuffer[floor(x)][floor(y)]){
	
				if(depthPoints.at(floor(i)) < 0.9){
					window.setPixelColour(floor(x), floor(y), colour);
					newDepthBuffer[floor(x)][floor(y)] = depthPoints.at(floor(i));
				}
				
			}
		}
		
		
	}
	return newDepthBuffer;
	
}



std::vector<CanvasPoint> interpolateCanvasPoints(CanvasPoint p1, CanvasPoint p2, int n){
	std::vector<CanvasPoint> result;
	//if(n==1)n=2;
	float xGap = (p2.x-p1.x)/(n-1);
	float yGap = (p2.y-p1.y)/(n-1);
	float zGap = (p2.depth-p1.depth)/(n-1);
	for (int i = 0; i < n; i++){
		float xNew = p1.x + xGap*i;
		float yNew = p1.y + yGap*i;
		float zNew = p1.depth + zGap*i;
		CanvasPoint newPoint;
		newPoint.x = xNew;
		newPoint.y = yNew;
		newPoint.depth = zNew;
		result.push_back(newPoint);
	}
	return result;
}

std::vector<TexturePoint> interpolateTexturePointsPerspective(CanvasPoint p1, CanvasPoint p2, int n){
	// get vertex with greatest depth
	float z0 = 0;
	float z1 = 0;

	float c0 = 0;
	float c1 = 0;

	float q = 0;

	float c = 0;

	p1.depth = 1;
	p2.depth = 1;

	std::vector<TexturePoint> result;
	TexturePoint t1 = p1.texturePoint;
	TexturePoint t2 = p2.texturePoint;
	if(n==1)n=2;
	if (n == 0) n = 2;
	
	float xGap = 0;

	xGap = (t2.x-t1.x)/(n-1);
	float yGap = (p2.texturePoint.y-p1.texturePoint.y)/(n-1); // may need to switch

	for (int i = 0; i < n; i++){

		q = i;

		// non linear interpolation of texture row
		if(p1.depth > p2.depth){ z0 = p1.depth; z1 = p2.depth; c0 = p1.texturePoint.y; c1 = p2.texturePoint.y;}
		else { z0 = p2.depth; z1 = p1.depth; c0 = p2.texturePoint.y ; c1 = p1.texturePoint.y; }

		float numerator = (c0/z0)*(1.0f-q) + (c1/z1)*q;
		float denominator = (1.0f/z0)*(1.0f-q) + (1.0f/z1)*q;

		c = numerator / denominator;

		std::cout << "c0: " << c0 <<  " c1: " << c1 << std::endl;

		TexturePoint newPoint;
		float xNew = t1.x + xGap*i;
		newPoint.x = xNew;
		newPoint.y = c;
		float yNew = t1.y + yGap*i;
		TexturePoint newPointNormal;
		newPointNormal.y = yNew;
		newPointNormal.x = xNew;

		if(c > 1000){
			std::cout << "LARGE C q: " << q <<  " z0: " << z0 << " z1: "  << z1 << std::endl;
		}

		std::cout << "perspective: " << newPoint << std::endl;
		std::cout << "normal: " << newPointNormal << std::endl;
		result.push_back(newPoint);

	}
	return result;


	

}

std::vector<TexturePoint> interpolateTexturePoints(TexturePoint t1, TexturePoint t2, int n){
	
	std::vector<TexturePoint> result;
	if(n==1)n=2;
	float xGap = (t2.x-t1.x)/(n-1);
	float yGap = (t2.y-t1.y)/(n-1);
//	std::cout << "t1: " << t1 << " t2: " << t2 << std::endl;
//	std::cout << "xGap: " << xGap << " yGap: " << yGap << std::endl;
	for (int i = 0; i < n; i++){
		float xNew = t1.x + xGap*i;
		float yNew = t1.y + yGap*i;
		TexturePoint newPoint;
		newPoint.x = xNew;
		newPoint.y = yNew;
		result.push_back(newPoint);
	//	std::cout << "normal: " << newPoint << std::endl;

	}
	

	return result;
}
std::vector<std::vector<uint32_t>> getPixelMatrixNormal(NormalMap &tMap){
	int width = tMap.width;
	int height = tMap.height;
	std::vector<std::vector<uint32_t>> pixelMatrix;

	// across x, down y

	// for y < height, define new vector
	for(int y = 0; y < height; y++){
		std::vector<uint32_t> rowValues;
		for(int x = 0; x < width; x++){
			
			// read x values from the map and add to rowValues
			rowValues.push_back(tMap.pixels.at(y*width + x));
		}
		pixelMatrix.push_back(rowValues);
	}
	// return pixelmatrix
	return pixelMatrix;
}

std::vector<std::vector<uint32_t>> getPixelMatrix(TextureMap &tMap){
	int width = tMap.width;
	int height = tMap.height;
	std::vector<std::vector<uint32_t>> pixelMatrix;

	// across x, down y

	// for y < height, define new vector
	for(int y = 0; y < height; y++){
		std::vector<uint32_t> rowValues;
		for(int x = 0; x < width; x++){
			
			// read x values from the map and add to rowValues
			rowValues.push_back(tMap.pixels.at(y*width + x));
		}
		pixelMatrix.push_back(rowValues);
	}
	// return pixelmatrix
	return pixelMatrix;
}
float maxXYdist(CanvasPoint &p1, CanvasPoint &p2){
	float stepCount=0;
	float xDiff = abs(p1.x - p2.x);
	float yDiff = abs(p1.y - p2.y);
	
	float xDiffSigned = p1.x - p2.x;
	float yDiffSigned = p1.y - p2.y;

	if(xDiff  > yDiff) stepCount = xDiff;
	else stepCount = yDiff;
	return stepCount;
}
std::vector<uint32_t> textureLineColours(CanvasPoint &p1, CanvasPoint &p2, TextureMap &tMap, std::vector<std::vector<uint32_t>> &pixelMatrix){


	

	std::vector<uint32_t> result;
	std::vector<TexturePoint> texturePoints;
	std::vector<CanvasPoint> canvasPoints;
	// texture point info

	// canvas point info
	float stepCountC = abs(p1.x- p2.x);
	if(stepCountC < 1)stepCountC = 1;
//	std::cout << "stepcountc " << stepCountC << std::endl;
	
	// get a list of the texture points p1 to p2. Get the acmount of pixels required by the line length

	// these are the coordinates of the texture we want

	texturePoints = interpolateTexturePoints(p2.texturePoint, p1.texturePoint, stepCountC);

	//now get the colours from the texture points
	


	for(int i = 0; i < floor(stepCountC)-1; i++){
		TexturePoint p = texturePoints.at(i);
		int x = floor(p.x);
		int y = floor(p.y);
		//std::cout << "accessing pixel matrix at "<< x << "," << y << std::endl;
		if(y>0 && x > 0){
			uint32_t colourToAdd = pixelMatrix[y][x];
			result.push_back(colourToAdd);
		}
		
	}
	std::cout << "got list of colours with size "<< result.size() << std::endl;
	return result;

}
// Canvas points require texture points
std::vector<std::vector<float>> drawLineTextureDepthBuffered(CanvasPoint &point1, CanvasPoint &point2, TextureMap &tMap ,DrawingWindow &window, std::vector<std::vector<uint32_t>> &pixelMatrix, std::vector<std::vector<float>> depthBuffer){

	std::vector<std::vector<float>> newDepthBuffer = depthBuffer;

	float xDiff = abs(point1.x - point2.x);
	float yDiff = abs(point1.y - point2.y);
	
	float xDiffSigned = point1.x - point2.x;
	float yDiffSigned = point1.y - point2.y;
	
	
	float stepCount= 0;
	if(xDiff  > yDiff) stepCount = xDiff;
	else stepCount = yDiff;
	if(stepCount <1) stepCount = 1;
	
	float xStepSize = xDiffSigned/stepCount;
	float yStepSize = yDiffSigned/stepCount;


	// draw a line from t1 to t2 on the texture map
//	std::cout << "begin textureLineColours" << std::endl;
//	std::cout << "----------------------" << std::endl;
//	std::cout << "point 1 " <<  point1.texturePoint <<  " point 2"  << point2.texturePoint << std::endl;
	std::vector<uint32_t> lineColours = textureLineColours(point1, point2, tMap, pixelMatrix);
	std::vector<float> depthPoints = interpolateDepthPoints(point1, point2, lineColours.size());
	std::cout << "done interpolateDepthPoints with size: "<< depthPoints.size() << std::endl;
	std::cout << "done lineColours with size: "<< lineColours.size() << std::endl;
//	std::cout << "stepCount is "<< stepCount << std::endl;
//	std::cout << "stepCount length" << stepCount <<std::endl;

	for (float i = 0.0; i<lineColours.size(); i++){
		float x = point2.x + (xStepSize*i);
		float y = point2.y + (yStepSize*i);

		/*
		float red = col.red;
		float green = col.green;
		float blue = col.blue;*/
//		std::cout << "final i is" << i <<std::endl;
	//	std::cout << "got here" << std::endl;
	//	std::cout << lineColours.size() << std::endl;
	//	std::cout << x << " , " << y <<std::endl;
		if(lineColours.size() > i && y<HEIGHT && x<WIDTH && x>=0 && y>=0){
			uint32_t colour = lineColours.at(floor(i));
			if(depthPoints.at(i) > depthBuffer[x][y]){
				window.setPixelColour(floor(x), floor(y), colour);
				depthBuffer[x][y] = depthPoints.at(i);
			}
			
		}

		
	}
	std::cout << "drawn line" << std::endl;
	return depthBuffer;
}

std::vector<std::vector<float>> drawLineTexture(CanvasPoint &point1, CanvasPoint &point2, TextureMap &tMap ,DrawingWindow &window, std::vector<std::vector<uint32_t>> &pixelMatrix){


	float xDiff = abs(point1.x - point2.x);
	float yDiff = abs(point1.y - point2.y);
	
	float xDiffSigned = point1.x - point2.x;
	float yDiffSigned = point1.y - point2.y;
	
	
	float stepCount= 0;
	if(xDiff  > yDiff) stepCount = xDiff;
	else stepCount = yDiff;
	if(stepCount <1) stepCount = 1;
	
	float xStepSize = xDiffSigned/stepCount;
	float yStepSize = yDiffSigned/stepCount;


	// draw a line from t1 to t2 on the texture map
	std::cout << "begin textureLineColours" << std::endl;
	std::cout << "----------------------" << std::endl;
	std::cout << "point 1 " <<  point1.texturePoint <<  " point 2"  << point2.texturePoint << std::endl;
	std::vector<uint32_t> lineColours = textureLineColours(point1, point2, tMap, pixelMatrix);
	std::vector<float> depthPoints = interpolateDepthPoints(point1, point2, stepCount);
	std::cout << "done textureLineColours with size: "<< lineColours.size() << std::endl;
	std::cout << "stepCount is "<< stepCount << std::endl;
//	std::cout << "stepCount length" << stepCount <<std::endl;

	for (float i = 0.0; i<lineColours.size(); i++){
		float x = point2.x + (xStepSize*i);
		float y = point2.y + (yStepSize*i);

		/*
		float red = col.red;
		float green = col.green;
		float blue = col.blue;*/
//		std::cout << "final i is" << i <<std::endl;
//		std::cout << "got here" << std::endl;
//		std::cout << lineColours.size() << std::endl;
		if(lineColours.size() > i){
			uint32_t colour = lineColours.at(floor(i));
			window.setPixelColour(floor(x), floor(y), colour);
			
		}

		
	}
}


CanvasPoint triangleMidIntersection(CanvasPoint yMax, CanvasPoint yMid, CanvasPoint yMin){
	//get y pos of y mid
	float interesectionY = yMid.y;
	
	//find x pos of the intersection with yMax, yMin
	// get inverse gradient of yMin yMax (dx/dy)
	float gradient = (yMin.x-yMax.x)/(yMin.y-yMax.y);
	float intersectionX = yMin.x + gradient*(interesectionY-yMin.y);



	//Create a new CanvasPoint at this intersection
	CanvasPoint flatBottomIntersection;
	flatBottomIntersection.x = intersectionX;
	flatBottomIntersection.y = interesectionY;


	// 
	float gradientZ = (yMin.depth-yMax.depth)/(yMin.y-yMax.y);
	float intersectionZ = yMin.depth + gradientZ*(interesectionY-yMin.y);



	//Create a new CanvasPoint at this intersection
//	CanvasPoint flatBottomIntersection;
	flatBottomIntersection.x = intersectionX;
	flatBottomIntersection.y = interesectionY;
	flatBottomIntersection.depth = intersectionZ;

	
	

	return flatBottomIntersection;
}


float max(float a, float b){
	if(a>b) return a;
	else return b;
}

void drawTexturedTriangle(CanvasTriangle &triangleVertices, DrawingWindow &window, TextureMap &tMap){

	//TextureMap tMap(texturePath);
	
	
	
	CanvasPoint p1 = triangleVertices.vertices[0];
	CanvasPoint p2 = triangleVertices.vertices[1];
	CanvasPoint p3 = triangleVertices.vertices[2];


	
	// now order the vertices from top to bottom.
	CanvasPoint yMax;
	CanvasPoint yMid;
	CanvasPoint yMin;
	
	// sort y pos order
	if(p1.y < p2.y && p1.y < p3.y){
		yMax = p1;
		if(p2.y<p3.y){yMid = p2; yMin = p3;}
		else {yMid = p3; yMin = p2;}
	} 
	else if(p2.y< p1.y && p2.y < p3.y){
		yMax = p2;
		if(p1.y<p3.y){yMid = p1; yMin = p3;}
		else {yMid = p3; yMin = p1;}
	} 
	else{
		yMax = p3;
		if(p1.y<p2.y){yMid = p1; yMin = p2;}
		else {yMid = p3; yMin = p2;}
	} 
//	std::cout << "Texturepoint yMax " <<yMax.texturePoint << std::endl;
//	std::cout << "Texturepoint yMid " <<yMid.texturePoint << std::endl;
//	std::cout << "Texturepoint yMin " <<yMin.texturePoint << std::endl;
	// now we need to map a texture point to each of these canvas points, including the interection
	
	// step 1 - get interection
	CanvasPoint yIntersection;
	yIntersection = triangleMidIntersection(yMax, yMid, yMin);
	// step 2 - assign texture point to the intersection
	// found by interpolating between yMax.texturePoint and yMin.texturePoint 
	
	// 
	float steps = abs(yMax.y-yMid.y)+2;
	std::vector<TexturePoint> sideTexturepointsR1 = interpolateTexturePoints(yMax.texturePoint, yMin.texturePoint, steps);
	
	// now get the nth element, equal to steps - (yMin.y - yIntersection.y)
	//int index = (int)round(steps - (abs(yMin.y - yIntersection.y)));
	float t = (abs(abs(yMax.y-yIntersection.y)))/(abs(yMax.y-yMin.y));
	float index = round(t*abs(yMax.y-yMin.y)); //---------------------------------
//	std::cout << "index" << index <<std::endl;
//	std::cout << "t" << t <<std::endl;
	TexturePoint intersectionTexturePoint = sideTexturepointsR1.at(index);
	
	yIntersection.texturePoint.x = intersectionTexturePoint.x;
	yIntersection.texturePoint.y = intersectionTexturePoint.y;


	
	// get the pixel matrix 
	std::vector<std::vector<uint32_t>> pixelMatrix = getPixelMatrix(tMap);
	
	// we have now got all of the relevant texturepoints for each vertex
	// draw the triangle like a raster
	
	int steps1 = abs(yMid.y-yMin.y)+2;
	std::vector<CanvasPoint> sideVectorsR1 = interpolateCanvasPoints(yMax, yMid, steps1);
	// gets list of texture points from yMax to yMid
	std::vector<TexturePoint> sideVectorsR1tex = interpolateTexturePoints(yMax.texturePoint, yMid.texturePoint, steps1);
	
	std::vector<CanvasPoint> sideVectorsL1 = interpolateCanvasPoints(yMax, yIntersection, steps1);
	// gets list of texture points from yMax to yIntersection
	std::vector<TexturePoint> sideVectorsL1tex = interpolateTexturePoints(yMax.texturePoint, yIntersection.texturePoint, steps1);
	
	for(float y = 0; y < steps1; y += 1.0){
		//draw line from edge to edge
		int i = (int)y;
		CanvasPoint point1 = sideVectorsR1.at(i);
		point1.texturePoint = sideVectorsR1tex.at(i);
		CanvasPoint point2 = sideVectorsL1.at(i);
		point2.texturePoint = sideVectorsL1tex.at(i);
		
		drawLineTexture(point1, point2, tMap, window, pixelMatrix);
		//draw a line using pixel colour interpolation
	}
	
	int steps2 = max(maxXYdist(yMid, yMin), maxXYdist(yIntersection, yMin));
	std::vector<CanvasPoint> sideVectorsR2 = interpolateCanvasPoints(yMid, yMin, steps2);
	std::vector<TexturePoint> sideVectorsR2tex = interpolateTexturePoints(yMid.texturePoint, yMin.texturePoint, steps2);
	
	std::vector<CanvasPoint> sideVectorsL2 = interpolateCanvasPoints(yIntersection, yMin, steps2);
	std::vector<TexturePoint> sideVectorsL2tex = interpolateTexturePoints(yIntersection.texturePoint, yMin.texturePoint, steps2);
	
	for(float y = 0; y < steps2; y += 1.0){
		//draw line from edge to edge
		int i = (int)y;
		CanvasPoint point1 = sideVectorsR2.at(i);
		point1.texturePoint = sideVectorsR2tex.at(i);
		CanvasPoint point2 = sideVectorsL2.at(i);
		point2.texturePoint = sideVectorsL2tex.at(i);
//		std::cout << "drawing line :" << point1 << point2 << std::endl;
		drawLineTexture(point1, point2, tMap, window, pixelMatrix);
		//draw a line using pixel colour interpolation
	}
	

	


 
}
std::vector<std::vector<float>> drawTexturedTriangleDepthBuffered(CanvasTriangle &triangleVertices, DrawingWindow &window, TextureMap &tMap, std::vector<std::vector<float>> depthBuffer){

	//TextureMap tMap(texturePath);
	std::vector<std::vector<float>> newDepthBuffer = depthBuffer;
	
	
	CanvasPoint p1 = triangleVertices.vertices[0];
	CanvasPoint p2 = triangleVertices.vertices[1];
	CanvasPoint p3 = triangleVertices.vertices[2];


	
	// now order the vertices from top to bottom.
	CanvasPoint yMax;
	CanvasPoint yMid;
	CanvasPoint yMin;
	
	// sort y pos order
	if(p1.y < p2.y && p1.y < p3.y){
		yMax = p1;
		if(p2.y<p3.y){yMid = p2; yMin = p3;}
		else {yMid = p3; yMin = p2;}
	} 
	else if(p2.y< p1.y && p2.y < p3.y){
		yMax = p2;
		if(p1.y<p3.y){yMid = p1; yMin = p3;}
		else {yMid = p3; yMin = p1;}
	} 
	else{
		yMax = p3;
		if(p1.y<p2.y){yMid = p1; yMin = p2;}
		else {yMid = p2; yMin = p1;}
	}
	std::cout << "Texturepoint yMax " <<yMax.texturePoint << std::endl;
	std::cout << "Texturepoint yMid " <<yMid.texturePoint << std::endl;
	std::cout << "Texturepoint yMin " <<yMin.texturePoint << std::endl;
	std::cout << "Depth yMax " <<yMax.depth << std::endl;
	std::cout << "Depth yMid " <<yMid.depth << std::endl;
	std::cout << "Depth yMin " <<yMin.depth << std::endl;
	// now we need to map a texture point to each of these canvas points, including the interection
	
	// step 1 - get interection
	CanvasPoint yIntersection;
	yIntersection = triangleMidIntersection(yMax, yMid, yMin);
	// step 2 - assign texture point to the intersection
	// found by interpolating between yMax.texturePoint and yMin.texturePoint 
	
	// 
	float steps = abs(yMax.y-yMid.y)+2;
	std::vector<TexturePoint> sideTexturepointsR1 = interpolateTexturePoints(yMax.texturePoint, yMin.texturePoint, steps);
	//std::vector<TexturePoint> sideTexturepointsR1 = interpolateTexturePointsPerspective(yMax, yMin, steps);
	
	// now get the nth element, equal to steps - (yMin.y - yIntersection.y)

	float t = (abs(abs(yMax.y-yIntersection.y)))/(abs(yMax.y-yMin.y));
	float index = round(t*abs(yMax.y-yMin.y)); //---------------------------------
	std::cout << "index" << index <<std::endl;
	std::cout << "t" << t <<std::endl;
	TexturePoint intersectionTexturePoint = sideTexturepointsR1.at(index);
	
	yIntersection.texturePoint.x = intersectionTexturePoint.x;
	yIntersection.texturePoint.y = intersectionTexturePoint.y;


	
	// get the pixel matrix 
	std::vector<std::vector<uint32_t>> pixelMatrix = getPixelMatrix(tMap);
	
	// we have now got all of the relevant texturepoints for each vertex
	// draw the triangle like a raster
	
	int steps1 = abs(yMax.y-yMid.y)+2;
	std::vector<CanvasPoint> sideVectorsR1 = interpolateCanvasPoints(yMax, yMid, steps1);
	// gets list of texture points from yMax to yMid
	std::vector<TexturePoint> sideVectorsR1tex = interpolateTexturePoints(yMax.texturePoint, yMid.texturePoint, steps1);
//	std::vector<TexturePoint> sideVectorsR1tex = interpolateTexturePointsPerspective(yMax, yMid, steps1);
	
	std::vector<CanvasPoint> sideVectorsL1 = interpolateCanvasPoints(yMax, yIntersection, steps1);
	// gets list of texture points from yMax to yIntersection
	std::vector<TexturePoint> sideVectorsL1tex = interpolateTexturePoints(yMax.texturePoint, yIntersection.texturePoint, steps1);
//	std::vector<TexturePoint> sideVectorsL1tex = interpolateTexturePointsPerspective(yMax, yIntersection, steps1);

	std::cout << "steps1 " << steps1 << std::endl;
	for(float y = 0; y < steps1; y += 1.0){
		//draw line from edge to edge
		int i = (int)y;
		CanvasPoint point1 = sideVectorsR1.at(i);
		point1.texturePoint = sideVectorsR1tex.at(i);
		CanvasPoint point2 = sideVectorsL1.at(i);
		point2.texturePoint = sideVectorsL1tex.at(i);
		std::cout << "drawLineTextureDepthBuffered " << point1 << point2 << std::endl;
		newDepthBuffer = drawLineTextureDepthBuffered(point1, point2, tMap, window, pixelMatrix, newDepthBuffer);
		//draw a line using pixel colour interpolation
	}
	
	int steps2 = abs(yMid.y-yMin.y)+2;
	std::vector<CanvasPoint> sideVectorsR2 = interpolateCanvasPoints(yMid, yMin, steps2);
	std::vector<TexturePoint> sideVectorsR2tex = interpolateTexturePoints(yMid.texturePoint, yMin.texturePoint, steps2);
	
	std::vector<CanvasPoint> sideVectorsL2 = interpolateCanvasPoints(yIntersection, yMin, steps2);
	std::vector<TexturePoint> sideVectorsL2tex = interpolateTexturePoints(yIntersection.texturePoint, yMin.texturePoint, steps2);
	
	for(float y = 0; y < steps2; y += 1.0){
		//draw line from edge to edge
		int i = (int)y;
		CanvasPoint point1 = sideVectorsR2.at(i);
		point1.texturePoint = sideVectorsR2tex.at(i);
		CanvasPoint point2 = sideVectorsL2.at(i);
		point2.texturePoint = sideVectorsL2tex.at(i);
		//		std::cout << "drawing line :" << point1 << point2 << std::endl;
		std::cout << "drawLineTextureDepthBuffered2 " << point1 << point2 << std::endl;
		newDepthBuffer = drawLineTextureDepthBuffered(point1, point2, tMap, window, pixelMatrix, newDepthBuffer);
		//draw a line using pixel colour interpolation
	}

	return newDepthBuffer;
	
}
std::vector<std::vector<float>> drawRasterisedTriangleDepthBuffered(CanvasTriangle &triangleVertices, Colour &col, DrawingWindow &window, std::vector<std::vector<float>> depthBuffer){
	std::vector<std::vector<float>> newDepthBuffer = depthBuffer;
	CanvasPoint p1 = triangleVertices.vertices[0];
	CanvasPoint p2 = triangleVertices.vertices[1];
	CanvasPoint p3 = triangleVertices.vertices[2];
	// now order the vertices from top to bottom.
	CanvasPoint yMax;
	CanvasPoint yMid;
	CanvasPoint yMin;
	if(p1.y < p2.y && p1.y < p3.y){
		yMax = p1;
		if(p2.y<p3.y){yMid = p2; yMin = p3;}
		else {yMid = p3; yMin = p2;}
	} 
	else if(p2.y< p1.y && p2.y < p3.y){
		yMax = p2;
		if(p1.y<p3.y){yMid = p1; yMin = p3;}
		else {yMid = p3; yMin = p1;}
	} 
	else{
		yMax = p3;
		if(p1.y<p2.y){yMid = p1; yMin = p2;}
		else {yMid = p2; yMin = p1;}
	}
	//Create a new CanvasPoint at this intersection

	CanvasPoint flatBottomIntersection;
	flatBottomIntersection = triangleMidIntersection(yMax, yMid, yMin);

	float lineMultiplier = 1;
	//Triangle 1 - yMax, yMid, flatBottomIntersection
	//int steps1 = max(maxXYdist(yMax, yMid), maxXYdist(yMax, flatBottomIntersection))*1.2;
	int steps1 = abs(yMax.y-yMid.y)+2;
	
	std::vector<CanvasPoint> sideVectorsR1 = interpolateCanvasPoints(yMax, yMid, steps1*lineMultiplier);
	std::vector<CanvasPoint> sideVectorsL1 = interpolateCanvasPoints(yMax, flatBottomIntersection, steps1*lineMultiplier);
	
	for(float y = 0; y < steps1*lineMultiplier; y += 1.0){
		//draw line from edge to edge
		int i = (int)y;
		newDepthBuffer = drawLineDepthBuffered(sideVectorsL1.at(i), sideVectorsR1.at(i), col, window, newDepthBuffer, lineMultiplier);
	}
	
	//Triangle 2 - flatBottomIntersection, yMid, yMin
	//int steps2 = max(maxXYdist(yMid, yMin), maxXYdist(flatBottomIntersection, yMin))*1.2;
	int steps2 = abs(yMid.y-yMin.y)+2;

	std::vector<CanvasPoint> sideVectorsR2;
	sideVectorsR2 = interpolateCanvasPoints(yMid, yMin, steps2*lineMultiplier);

	std::vector<CanvasPoint> sideVectorsL2;
	sideVectorsL2 = interpolateCanvasPoints(flatBottomIntersection, yMin, steps2*lineMultiplier);

	for(float y = 0; y < steps2*lineMultiplier; y += 1.0){
		//draw line from edge to edge
		int i = (int)y;
		newDepthBuffer = drawLineDepthBuffered(sideVectorsL2.at(i), sideVectorsR2.at(i), col, window, newDepthBuffer, lineMultiplier);
	}
	return newDepthBuffer;
}

void drawRasterisedTriangle(CanvasTriangle &triangleVertices, Colour &col, DrawingWindow &window){
	CanvasPoint p1 = triangleVertices.vertices[0];
	CanvasPoint p2 = triangleVertices.vertices[1];
	CanvasPoint p3 = triangleVertices.vertices[2];
	// now order the vertices from top to bottom.
	CanvasPoint yMax;
	CanvasPoint yMid;
	CanvasPoint yMin;
	if(p1.y < p2.y && p1.y < p3.y){
		yMax = p1;
		if(p2.y<p3.y){yMid = p2; yMin = p3;}
		else {yMid = p3; yMin = p2;}
	} 
	else if(p2.y< p1.y && p2.y < p3.y){
		yMax = p2;
		if(p1.y<p3.y){yMid = p1; yMin = p3;}
		else {yMid = p3; yMin = p1;}
	} 
	else{
		yMax = p3;
		if(p1.y<p2.y){yMid = p1; yMin = p2;}
		else {yMid = p2; yMin = p1;}
	}
	//Create a new CanvasPoint at this intersection
	std::cout << "points:" << std::endl;
	std::cout << yMax << yMid << yMin << std::endl;
	CanvasPoint flatBottomIntersection;
	flatBottomIntersection = triangleMidIntersection(yMax, yMid, yMin);

	//Triangle 1 - yMax, yMid, flatBottomIntersection
	int steps1 = max(maxXYdist(yMax, yMid), maxXYdist(yMax, flatBottomIntersection));
	
	std::vector<CanvasPoint> sideVectorsR1 = interpolateCanvasPoints(yMax, yMid, steps1*1.1);
	std::vector<CanvasPoint> sideVectorsL1 = interpolateCanvasPoints(yMax, flatBottomIntersection, steps1*1.1);
	
	for(float y = 0; y < steps1*1.1 - 1; y += 1.0){
		//draw line from edge to edge
		int i = (int)y;
		drawLine(sideVectorsL1.at(i), sideVectorsR1.at(i), col, window);
	}
	
	//Triangle 2 - flatBottomIntersection, yMid, yMin
	int steps2 = max(maxXYdist(yMid, yMin), maxXYdist(flatBottomIntersection, yMin));
	std::cout << "steps2:" << steps2 <<  std::endl;
	std::vector<CanvasPoint> sideVectorsR2;
	sideVectorsR2 = interpolateCanvasPoints(yMid, yMin, steps2*1.1);

	std::vector<CanvasPoint> sideVectorsL2;
	sideVectorsL2 = interpolateCanvasPoints(flatBottomIntersection, yMin, steps2*1.1);

	for(float y = 0; y < steps2*1.1 -1; y += 1.0){
		//draw line from edge to edge
		int i = (int)y;
		drawLine(sideVectorsL2.at(i), sideVectorsR2.at(i), col, window);
	}
}
// requires depth points to be set before calling
std::vector<std::vector<float>> drawDepthBufferedStrokedTriangle(CanvasTriangle &triangleVertices, Colour &col, DrawingWindow &window, std::vector<std::vector<float>> depthBuffer){
	std::vector<std::vector<float>> newDepthBuffer = depthBuffer;
	
	CanvasPoint p1 = triangleVertices.vertices[0];
	CanvasPoint p2 = triangleVertices.vertices[1];
	CanvasPoint p3 = triangleVertices.vertices[2];

	newDepthBuffer = drawLineDepthBuffered(p1, p2, col, window, newDepthBuffer, 1);
	newDepthBuffer = drawLineDepthBuffered(p2, p3, col, window, newDepthBuffer, 1);
	newDepthBuffer = drawLineDepthBuffered(p1, p3, col, window, newDepthBuffer, 1);

	return newDepthBuffer;
}



void drawStrokedTriangle(CanvasTriangle &triangleVertices, Colour &col, DrawingWindow &window){
	CanvasPoint p1 = triangleVertices.vertices[0];
	CanvasPoint p2 = triangleVertices.vertices[1];
	CanvasPoint p3 = triangleVertices.vertices[2];
	// draw line from 1 to 2
	drawLine(p1, p2, col, window);
	// draw line from 2 to 3
	drawLine(p2, p3, col, window);
	// draw line from 1 to 3
	drawLine(p1, p3, col, window);
}
/*
void draw(DrawingWindow &window) {
	//window.clearPixels();
	/*for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = rand() % 256;
			float green = 0.0;
			float blue = 0.0;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}

	// Drawing lines test


}*/


/*
void update(DrawingWindow &window) {
	// Function for performing animation (shifting artifacts or moving the camera)
}*/

void generateAndDrawTriangle(DrawingWindow &window, int typeID){
	CanvasTriangle trianglePoints;
	
	CanvasPoint p1;
	CanvasPoint p2;
	CanvasPoint p3;
	p1.x = float(rand() % window.width);
	p1.y = float(rand() % window.height);
	

	p2.x = float(rand() % window.width);
	p2.y = float(rand() % window.height);

	p3.x = float(rand() % window.width);
	p3.y = float(rand() % window.height);
	
	trianglePoints.vertices[0] = p1;
	trianglePoints.vertices[1] = p2;
	trianglePoints.vertices[2] = p3;
//	std::cout << trianglePoints << std::endl;
	Colour col;
	col.red = rand()%256;
	col.green = rand()%256;
	col.blue = rand()%256;

	if(typeID == 0) drawStrokedTriangle(trianglePoints, col, window);
	else if (typeID == 1) drawRasterisedTriangle(trianglePoints, col, window);
	//else if (typeID == 2) drawTexturedTriangle(trianglePoints, window);
}
void fixedTriangle(DrawingWindow &window){
	CanvasTriangle trianglePoints;
	CanvasTriangle texturePoints;
	CanvasPoint p1;
	CanvasPoint p2;
	CanvasPoint p3;

	TexturePoint t1;
	TexturePoint t2;
	TexturePoint t3;

	p1.x = 160;
	p1.y = 10;
	
	p2.x = 300;
	p2.y = 230;

	p3.x = 10;
	p3.y = 150;

	t1.x = 195;
	t1.y = 5;
	
	t2.x = 395;
	t2.y = 380;

	t3.x = 65;
	t3.y = 330;

	p1.texturePoint = t1;
	p2.texturePoint = t2;
	p3.texturePoint = t3;
	
	trianglePoints.vertices[0] = p1;
	trianglePoints.vertices[1] = p2;
	trianglePoints.vertices[2] = p3;

//	drawTexturedTriangle(trianglePoints, window);
}

/*void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
		else if (event.key.keysym.sym == SDLK_u) generateAndDrawTriangle(window, 0);
		else if (event.key.keysym.sym == SDLK_r) generateAndDrawTriangle(window, 1);
		else if (event.key.keysym.sym == SDLK_t) fixedTriangle(window);
	} else if (event.type == SDL_MOUSEBUTTONDOWN) window.savePPM("output.ppm");
}*/

/*int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		update(window);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}*/

