#include "Mesh.h"

Mesh::Mesh()
{
}

Mesh::~Mesh()
{
	vertices.clear();
	normals.clear();
	indices.clear();
	colors.clear();
}

void Mesh::scale(float scale)
{
	for(unsigned int i = 0; i < vertices.size(); ++i)
	{
		vertices[i] = vertices[i] * scale;
	}
}
void Mesh::setColor(glm::vec3 color)
{
	colors.resize(vertices.size());
	for(unsigned int i = 0; i < colors.size(); ++i)
	{
		colors[i] = color;
	}
}
int Mesh::readOBJ(const char* file)
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
	std::string index_string;
	while(!readFile.eof())
	{
		getline(readFile, line);
		std::istringstream iss(line);
		iss>>type;
		if(type[0] == '#')
			continue;
		if(!type.compare("v"))
		{
			iss>>vertex.x>>vertex.y>>vertex.z;
			vertices.push_back(vertex);
		}
		else if(!type.compare("f"))
		{
			iss>>index_string;
			index.x = atoi(index_string.substr(0, index_string.find("/", 0)).c_str()) - 1;
			iss>>index_string;
			index.y = atoi(index_string.substr(0, index_string.find("/", 0)).c_str()) - 1;
			iss>>index_string;
			index.z = atoi(index_string.substr(0, index_string.find("/", 0)).c_str()) - 1;
			indices.push_back(index.x);
			indices.push_back(index.y);
			indices.push_back(index.z);
		}
		else if(!type.compare("vn"))
		{
			iss>>normal.x>>normal.y>>normal.z;
			normals.push_back(normal);
		}
	}
	return 1;
}

int Mesh::writeOBJ(const char* file)
{
	return 1;
}

void Mesh::draw(const VBO& vbos)
{
	// set color
	setColor(glm::vec3(0.4, 0.9, 0.8));

	// vertices
	glBindBuffer(GL_ARRAY_BUFFER, vbos.m_vbo);
	glBufferData(GL_ARRAY_BUFFER, 3 * vertices.size() * sizeof(float), &vertices[0], GL_STREAM_DRAW);

	// colors
	glBindBuffer(GL_ARRAY_BUFFER, vbos.m_cbo);
	glBufferData(GL_ARRAY_BUFFER, 3 * colors.size() * sizeof(float), &colors[0], GL_STREAM_DRAW);

	// normals
	glBindBuffer(GL_ARRAY_BUFFER, vbos.m_nbo);
	glBufferData(GL_ARRAY_BUFFER, 3 * normals.size() * sizeof(float), &normals[0], GL_STREAM_DRAW);

	// indices
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos.m_ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned short), &indices[0], GL_STATIC_DRAW);

	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	glEnableVertexAttribArray(2);

	glBindBuffer(GL_ARRAY_BUFFER, vbos.m_vbo);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, vbos.m_cbo);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, vbos.m_nbo);
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos.m_ibo);
	glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_SHORT, 0);//GL_UNSIGNED_INT

	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
	glDisableVertexAttribArray(2);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}