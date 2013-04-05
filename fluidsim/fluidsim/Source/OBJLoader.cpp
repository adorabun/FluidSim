#include "OBJLoader.h"

OBJLoader::OBJLoader()
{
}

OBJLoader::~OBJLoader()
{
	vertices.clear();
	normals.clear();
	indices.clear();
}

int OBJLoader::readOBJ(const char* file)
{
	std::ifstream readFile;
	readFile.open(file);
	if(!readFile.is_open())
	{
		std::cerr<<"Cannot open OBJ file."<<std::endl;
		readFile.close();
		return 0;
	}
	std::string line;
	std::string type;
	glm::vec3 vertex, normal, index;
	while(!readFile.eof())
	{
		getline(readFile, line);
		std::istringstream iss(line);
		iss>>type;
		if(type[0] == '#')
			continue;
		if(type.compare("v"))
		{
			iss>>vertex.x>>vertex.y>>vertex.z;
			vertices.push_back(vertex);
		}
		else if(type.compare("f"))
		{
			iss>>index.x>>index.y>>index.z;
			indices.push_back(index.x);
			indices.push_back(index.y);
			indices.push_back(index.z);
		}
		else if(type.compare("vn"))
		{
			iss>>normal.x>>normal.y>>normal.z;
			normals.push_back(normal);
		}
	}
	return 1;
}

int OBJLoader::writeOBJ(const char* file)
{
	std::ofstream writeFile;

	return 1;
}