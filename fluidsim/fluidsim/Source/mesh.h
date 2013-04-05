#ifndef MESH_H_INCLUDED
#define MESH_H_INCLUDED

#include <vector>
#include "math_headers.h"
#include "openGL_headers.h"
#include <stdio.h>

class mesh{
public:
	mesh();
	void loadObj(const char* filename);
	void saveObj(const char* filename);
	void Draw(const VBO& vbos);
public:	
	std::vector<glm::vec4> vertices;
	std::vector<glm::vec3> normals;
	std::vector<unsigned short> indices;


}

#endif