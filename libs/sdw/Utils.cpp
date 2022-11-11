#include <algorithm>
#include <sstream>
#include "Utils.h"


std::vector<std::string> split(const std::string &line, char delimiter) {
	auto haystack = line;
	std::vector<std::string> tokens;
	size_t pos;
	while ((pos = haystack.find(delimiter)) != std::string::npos) {
		tokens.push_back(haystack.substr(0, pos));
		haystack.erase(0, pos + 1);
	}
	// Push the remaining chars onto the vector
	tokens.push_back(haystack);
	return tokens;
}

std::vector<glm::vec3> generatePointsInSphere(glm::vec3 origin, float radius, int n){
	std::vector<glm::vec3> result;


	
	int i = 6;
	// make i lights per layer
	// calculate how many layers in n
	float layerCount = (float) n/i;
	layerCount = ceil(layerCount);

	float curLayer = 1;

	int directionId = 0;

	for(int count = 0; count < n; count ++){
		// forward, back, left, right * radius * (curLayer/layerCount)
		glm::vec3 pos;
		if(directionId == 0) pos = origin + glm::vec3(0,0,1)*radius*(curLayer/layerCount);
		if(directionId == 1) pos = origin + glm::vec3(0,1,0)*radius*(curLayer/layerCount);
		if(directionId == 2) pos = origin + glm::vec3(1,0,0)*radius*(curLayer/layerCount);
		if(directionId == 3) pos = origin + glm::vec3(0,0,-1)*radius*(curLayer/layerCount);
		if(directionId == 4) pos = origin + glm::vec3(0,-1,0)*radius*(curLayer/layerCount);
		if(directionId == 5) pos = origin + glm::vec3(-1,0,0)*radius*(curLayer/layerCount);
		directionId ++;
		if(directionId > 5) directionId = 0;
		if(count >= i * curLayer){
			curLayer += 1;
		}
		result.push_back(pos);
	}
	

	return result;



}


