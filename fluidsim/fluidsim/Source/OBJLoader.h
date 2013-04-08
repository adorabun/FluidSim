#ifndef OBJ_LOADER_H_INCLUDED
#define OBJ_LOADER_H_INCLUDED

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include "math_headers.h"
#include "openGL_headers.h"

class OBJLoader
{
public:
	OBJLoader();
	~OBJLoader();
	int readOBJ(const char* file);
	int writeOBJ(const char* file);
	void Draw(const VBO& vbos);
private:
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> normals;
	std::vector<unsigned short> indices;
};

#endif