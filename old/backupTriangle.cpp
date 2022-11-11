#include <CanvasTriangle.h>
#include <TextureMap.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <CanvasPoint.h>
#include <Colour.h>

#define WIDTH 320
#define HEIGHT 240

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
		window.setPixelColour(int(x), int(y), colour);
	}
}

std::vector<CanvasPoint> interpolateCanvasPoints(CanvasPoint p1, CanvasPoint p2, int n){
	std::vector<CanvasPoint> result;
	if(n==1)n=2;
	float xGap = (p2.x-p1.x)/(n-1);
	float yGap = (p2.y-p1.y)/(n-1);
	for (int i = 0; i < n; i++){
		float xNew = p1.x + xGap*i;
		float yNew = p1.y + yGap*i;
		CanvasPoint newPoint;
		newPoint.x = xNew;
		newPoint.y = yNew;
		result.push_back(newPoint);
	}
	return result;
}

std::vector<TexturePoint> interpolateTexturePoints(TexturePoint t1, TexturePoint t2, int n){
	std::cout << "interpolateTexturePoints tp1 " << t1 << " tp2 " << t2 << std::endl;
	std::cout << "interpolateTexturePoints n " << n << std::endl;
	std::vector<TexturePoint> result;
	if(n==1)n=2;
	float xGap = (t2.x-t1.x)/(n-1);
	float yGap = (t2.y-t1.y)/(n-1);
	for (int i = 0; i < n; i++){
		float xNew = t1.x + xGap*i;
		float yNew = t1.y + yGap*i;
		TexturePoint newPoint;
		newPoint.x = xNew;
		newPoint.y = yNew;
		result.push_back(newPoint);
		std::cout << "interpolateTexturePoints tp is " << newPoint << std::endl;
	}
	std::cout << "got interpolateTexturePoints" << std::endl;
	
	return result;
}

std::vector<std::vector<uint32_t>> getPixelMatrix(TextureMap tMap){
	int width = tMap.width;
	int height = tMap.height;
	std::vector<std::vector<uint32_t>> pixelMatrix;

	// across x, down y

	// for y < height, define new vector
	for(int y = 0; y < height; y++){
		std::vector<uint32_t> rowValues;
		for(int x = 0; x < width; x++){
			
			// read x values from the map and add to rowValues
			rowValues.push_back(tMap.pixels.at(y*height + x));
		}
		pixelMatrix.push_back(rowValues);
	}
	// return pixelmatrix
	return pixelMatrix;
}
std::vector<uint32_t> textureLineColours(CanvasPoint &p1, CanvasPoint &p2, TextureMap &tMap, std::vector<std::vector<uint32_t>> &pixelMatrix){

//	std::cout << "p1 texture point:" << p1.texturePoint <<std::endl;
	

	std::vector<uint32_t> result;
	std::vector<TexturePoint> texturePoints;
	std::vector<CanvasPoint> canvasPoints;
	// texture point info

	// canvas point info
	float xDiffC = abs(p1.x - p2.x);
	float yDiffC = abs(p1.y - p2.y);
	
	float xDiffSignedC = p1.x - p2.x;
	float yDiffSignedC = p1.y - p2.y;
	float stepCountC= 0;
	if(xDiffC  > yDiffC) stepCountC = xDiffC;
	else stepCountC = yDiffC;

	
	// get a list of the texture points p1 to p2. Get the acmount of pixels required by the line length

	// these are the coordinates of the texture we want
	//std::cout << "begin interpolateTexturePoints:"  <<std::endl;
	texturePoints = interpolateTexturePoints(p1.texturePoint, p2.texturePoint, stepCountC);
	//std::cout << "begin interpolateCanvasPoints:"  <<std::endl;
	canvasPoints = interpolateCanvasPoints(p1, p2, stepCountC);
	//std::cout << "interpolateCanvasPoints success"  <<std::endl;

	//std::cout << "stepCountC"<< stepCountC  <<std::endl;
	//now get the colours from the texture points
	


	for(int i = 0; i < round(stepCountC)-1; i++){
		TexturePoint p = texturePoints.at(i);
		int x = round(p.x);
		int y = round(p.y);
		uint32_t colourToAdd = pixelMatrix[y][x];
		result.push_back(colourToAdd);
	}
	return result;

}
// Canvas points require texture points
void drawLineTexture(CanvasPoint &point1, CanvasPoint &point2, TextureMap &tMap ,DrawingWindow &window, std::vector<std::vector<uint32_t>> &pixelMatrix){

	
	//std::cout << pixelMatrix[0][0] << std::endl;
	float xDiff = abs(point1.x - point2.x);
	float yDiff = abs(point1.y - point2.y);
	
	float xDiffSigned = point1.x - point2.x;
	float yDiffSigned = point1.y - point2.y;
	
	
	float stepCount= 0;
	if(xDiff  > yDiff) stepCount = xDiff;
	else stepCount = yDiff;
	float xStepSize = xDiffSigned/stepCount;
	float yStepSize = yDiffSigned/stepCount;


	// draw a line from t1 to t2 on the texture map
	std::cout << "getting line colours" << std::endl;
	std::cout << "points" << point1 << point2 << std::endl;
	std::vector<uint32_t> lineColours = textureLineColours(point1, point2, tMap, pixelMatrix);
	std::cout << "got line colours" << std::endl;
//	std::cout << "stepCount length" << stepCount <<std::endl;
	
	for (float i = 0.0; i<stepCount-2; i++){
		float x = point2.x + (xStepSize*i);
		float y = point2.y + (yStepSize*i);
		std::cout << "got here" << lineColours.size() <<std::endl;
		/*
		float red = col.red;
		float green = col.green;
		float blue = col.blue;*/
//		std::cout << "final i is" << i <<std::endl;
		uint32_t colour = lineColours.at(round(i));
		window.setPixelColour(int(round(x)), int(round(y)), colour);
		
	}
}
CanvasPoint triangleMidIntersection(CanvasPoint yMax, CanvasPoint yMid, CanvasPoint yMin){
	//get y pos of y mid
	float interesectionY = yMid.y;
	//find x pos of the intersection with yMax, yMin
	// get inverse gradient of yMin yMax (dx/dy)
	float gradient = (yMax.x-yMin.x)/(yMax.y-yMin.y);
	float intersectionX = yMin.x + gradient*(interesectionY-yMin.y);

	//Create a new CanvasPoint at this intersection
	CanvasPoint flatBottomIntersection;
	flatBottomIntersection.x = intersectionX;
	flatBottomIntersection.y = interesectionY;
	return flatBottomIntersection;
}

void drawTexturedTriangle(CanvasTriangle &triangleVertices, CanvasTriangle &textureVertices, DrawingWindow &window){
	std::cout << "begin create tmap" << std::endl;
	TextureMap tMap("src/texture.ppm");
	
	
	
	CanvasPoint p1 = triangleVertices.vertices[0];
	CanvasPoint p2 = triangleVertices.vertices[1];
	CanvasPoint p3 = triangleVertices.vertices[2];

	CanvasPoint t1 = textureVertices.vertices[0];
	CanvasPoint t2 = textureVertices.vertices[1];
	CanvasPoint t3 = textureVertices.vertices[2];

	// map texture points to each triangle vertex(excluding intersection)
	p1.texturePoint.x = t1.x; p1.texturePoint.y = t1.y;
	p2.texturePoint.x = t2.x; p2.texturePoint.y = t2.y;
	p3.texturePoint.x = t3.x; p3.texturePoint.y = t3.y;
	
	// now order the vertices from top to bottom.
	CanvasPoint yMax;
	CanvasPoint yMid;
	CanvasPoint yMin;
	
	// sort y pos order
	if(p1.y > p2.y && p1.y > p3.y){
		yMax = p1;
		if(p2.y>p3.y){yMid = p2; yMin = p3;}
		else {yMid = p3; yMin = p2;}
	} 
	else if(p2.y> p1.y && p2.y > p3.y){
		yMax = p2;
		if(p1.y>p3.y){yMid = p1; yMin = p3;}
		else {yMid = p3; yMin = p1;}
	} 
	else{
		yMax = p3;
		if(p1.y>p2.y){yMid = p1; yMin = p2;}
		else {yMid = p3; yMin = p2;}
	} 
	// now we need to map a texture point to each of these canvas points, including the interection

	// step 1 - get interection
	CanvasPoint yIntersection;
	yIntersection = triangleMidIntersection(yMax, yMid, yMin);
	// step 2 - assign texture point to the intersection
	// found by interpolating between yMax.texturePoint and yMin.texturePoint 

	// 
	float steps = abs(yMax.y-yMin.y);
	std::vector<TexturePoint> sideTexturepointsR1 = interpolateTexturePoints(yMax.texturePoint, yMin.texturePoint, steps);
	
	// now get the nth element, equal to steps - (yMin.y - yIntersection.y)
	int index = (int)round(steps - (abs(yMin.y - yIntersection.y)));


	TexturePoint intersectionTexturePoint = sideTexturepointsR1.at(index);
	yIntersection.texturePoint.x = intersectionTexturePoint.x;
	yIntersection.texturePoint.y = intersectionTexturePoint.y;

	// get the pixel matrix 
	std::vector<std::vector<uint32_t>> pixelMatrix = getPixelMatrix(tMap);

	// we have now got all of the relevant texturepoints for each vertex
	// draw the triangle like a raster
	
	int steps1 = abs(yMax.y-yMid.y) + 1;
	std::vector<CanvasPoint> sideVectorsR1 = interpolateCanvasPoints(yMax, yMid, steps1);
	std::vector<TexturePoint> sideVectorsR1tex = interpolateTexturePoints(yMax.texturePoint, yMid.texturePoint, steps1);
	
	std::vector<CanvasPoint> sideVectorsL1 = interpolateCanvasPoints(yMax, yIntersection, steps1);
	std::vector<TexturePoint> sideVectorsL1tex = interpolateTexturePoints(yMax.texturePoint, yIntersection.texturePoint, steps1);

	for(float y = 0; y < steps1; y += 1.0){
		//draw line from edge to edge
		int i = (int)y;
		CanvasPoint point1 = sideVectorsR1.at(i);
		point1.texturePoint = sideVectorsR1tex.at(i);
		CanvasPoint point2 = sideVectorsL1.at(i);
		point2.texturePoint = sideVectorsL1tex.at(i);
//		std::cout << "drawing line :" << point1 << point2 << std::endl;
		drawLineTexture(point1, point2, tMap, window, pixelMatrix);
		//draw a line using pixel colour interpolation
	}
	int steps2 = abs(yMid.y-yMin.y) + 1;
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



void drawRasterisedTriangle(CanvasTriangle &triangleVertices, Colour &col, DrawingWindow &window){
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


	//Create a new CanvasPoint at this intersection
	
	CanvasPoint flatBottomIntersection;
	flatBottomIntersection = triangleMidIntersection(yMax, yMid, yMin);

	//Triangle 1 - yMax, yMid, flatBottomIntersection
	int steps1 = abs(yMax.y-yMid.y) + 1;
	
	std::vector<CanvasPoint> sideVectorsR1 = interpolateCanvasPoints(yMax, yMid, steps1);
	std::vector<CanvasPoint> sideVectorsL1 = interpolateCanvasPoints(yMax, flatBottomIntersection, steps1);
	
	for(float y = 0; y < steps1; y += 1.0){
		//draw line from edge to edge
		int i = (int)y;
		drawLine(sideVectorsL1.at(i), sideVectorsR1.at(i), col, window);
	}
	
	//Triangle 2 - flatBottomIntersection, yMid, yMin
	int steps2 = (yMid.y-yMin.y) + 1;
	std::vector<CanvasPoint> sideVectorsR2;
	sideVectorsR2 = interpolateCanvasPoints(yMid, yMin, steps2);

	std::vector<CanvasPoint> sideVectorsL2;
	sideVectorsL2 = interpolateCanvasPoints(flatBottomIntersection, yMin, steps2);

	for(float y = 0; y < steps2; y += 1.0){
		//draw line from edge to edge
		int i = (int)y;
		CanvasPoint point1 = sideVectorsL2.at(i);
		CanvasPoint point2 = sideVectorsR2.at(i);

		drawLine(point1, point2, col, window);
	}


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
	}*/

	// Drawing lines test

	
	float x1 = 100;
	float y1 = 100;
	CanvasPoint point1;
	point1.x = x1;
	point1.y = y1;

	float x2 = 200;
	float y2 = 200;
	CanvasPoint point2;
	point2.x = x2;
	point2.y = y2;
	Colour col;
	col.red = 255;
	drawLine(point1, point2, col ,window); 
}



void update(DrawingWindow &window) {
	// Function for performing animation (shifting artifacts or moving the camera)
}

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
	else if (typeID == 2) drawTexturedTriangle(trianglePoints, trianglePoints, window);
}
void fixedTriangle(DrawingWindow &window){
	CanvasTriangle trianglePoints;
	CanvasTriangle texturePoints;
	CanvasPoint p1;
	CanvasPoint p2;
	CanvasPoint p3;

	CanvasPoint t1;
	CanvasPoint t2;
	CanvasPoint t3;

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
	
	trianglePoints.vertices[0] = p1;
	trianglePoints.vertices[1] = p2;
	trianglePoints.vertices[2] = p3;
	texturePoints.vertices[0] = t1;
	texturePoints.vertices[1] = t2;
	texturePoints.vertices[2] = t3;
	drawTexturedTriangle(trianglePoints, texturePoints, window);
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
		else if (event.key.keysym.sym == SDLK_u) generateAndDrawTriangle(window, 0);
		else if (event.key.keysym.sym == SDLK_r) generateAndDrawTriangle(window, 1);
		else if (event.key.keysym.sym == SDLK_t) fixedTriangle(window);
	} else if (event.type == SDL_MOUSEBUTTONDOWN) window.savePPM("output.ppm");
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		update(window);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
