#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include "SDL.h"

#include "CanvasTriangle.h"
#include "TextureMap.h"
#include "NormalMap.h"
#include "DrawingWindow.h"
#include "Utils.h"
#include "fstream"
#include "vector"
#include "CanvasPoint.h"
#include "Colour.h"



	void drawLine(CanvasPoint &point1, CanvasPoint &point2, Colour &col ,DrawingWindow &window);
	void drawStrokedTriangle(CanvasTriangle &triangleVertices, Colour &col, DrawingWindow &window);
	void drawRasterisedTriangle(CanvasTriangle &triangleVertices, Colour &col, DrawingWindow &window);
	void drawTexturedTriangle(CanvasTriangle &triangleVertices, DrawingWindow &window, TextureMap &texturePath);
	void drawPixel(CanvasPoint &point, Colour &col, DrawingWindow &window);
	std::vector<std::vector<float>> drawDepthBufferedStrokedTriangle(CanvasTriangle &triangleVertices, Colour &col, DrawingWindow &window, std::vector<std::vector<float>> depthBuffer);
	std::vector<std::vector<float>> drawRasterisedTriangleDepthBuffered(CanvasTriangle &triangleVertices, Colour &col, DrawingWindow &window, std::vector<std::vector<float>> depthBuffer);
	std::vector<std::vector<float>> drawTexturedTriangleDepthBuffered(CanvasTriangle &triangleVertices, DrawingWindow &window, TextureMap &tMap, std::vector<std::vector<float>> depthBuffer);
	std::vector<std::vector<uint32_t>> getPixelMatrix(TextureMap &tMap);
	std::vector<std::vector<uint32_t>> getPixelMatrixNormal(NormalMap &tMap);
