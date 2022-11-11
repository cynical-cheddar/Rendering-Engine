#pragma once

#include <iostream>
#include <fstream>
#include <stdexcept>
#include "Utils.h"

class NormalMap {
public:
	size_t width;
	size_t height;
	std::vector<uint32_t> pixels;

	NormalMap();
	NormalMap(const std::string &filename);
	friend std::ostream &operator<<(std::ostream &os, const NormalMap &point);
};
