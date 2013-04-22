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
	face_normals.clear();
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
	glm::vec3 vertex, normal, index, normal_index, face_normal;
	std::string index_string;
	unsigned found;
	std::vector<glm::vec3> temp_normals;
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
		else if(!type.compare("vn"))
		{
			iss>>normal.x>>normal.y>>normal.z;
			temp_normals.push_back(normal);
			normals.push_back(normal);
		}
		else if(!type.compare("f"))
		{
			iss>>index_string;
			found = index_string.find("/", 0);
			index.x = atoi(index_string.substr(0, found).c_str()) - 1;
			found = index_string.find("/", found + 1);
			normal_index.x = atoi(index_string.substr(found + 1).c_str()) - 1;
			iss>>index_string;
			found = index_string.find("/", 0);
			index.y = atoi(index_string.substr(0, found).c_str()) - 1;
			found = index_string.find("/", found + 1);
			normal_index.y = atoi(index_string.substr(found + 1).c_str()) - 1;
			iss>>index_string;
			found = index_string.find("/", 0);
			index.z = atoi(index_string.substr(0, found).c_str()) - 1;
			found = index_string.find("/", found + 1);
			normal_index.z = atoi(index_string.substr(found + 1).c_str()) - 1;
			indices.push_back(index.x);
			indices.push_back(index.y);
			indices.push_back(index.z);
			// set vertex normal according to normal index
			normals[index.x] = temp_normals[normal_index.x];
			normals[index.y] = temp_normals[normal_index.y];
			normals[index.z] = temp_normals[normal_index.z];
			// calculate face normals
			if(vertices[index.x] == vertices[index.y] || vertices[index.x] == vertices[index.z] || vertices[index.y] == vertices[index.z])
			{
				face_normal = glm::normalize((normals[index.x] + normals[index.y] + normals[index.z]) / 3.0f);
				face_normals.push_back(face_normal);
			}
			else
			{
				face_normal = glm::normalize(glm::cross(vertices[index.y] - vertices[index.x], vertices[index.z] - vertices[index.x]));
				//face_normal = glm::normalize((normals[index.x] + normals[index.y] + normals[index.z]) / 3.0f);
				face_normals.push_back(face_normal);
			}
		}
	}
	temp_normals.clear();
	return 1;
}

int Mesh::writeOBJ(const char* file)
{
	return 1;
}

void Mesh::draw(const VBO& vbos)
{
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	// set color
	setColor(glm::vec3(0.6, 0.9, 0.8));

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

double triangleArea(glm::vec3 const& P1, glm::vec3 const& P2, glm::vec3 const& P3)
{
	double det1 = glm::determinant(glm::mat3(glm::vec3(P1.y, P1.z, 1), glm::vec3(P2.y, P2.z, 1), glm::vec3(P3.y, P3.z, 1)));
	double det2 = glm::determinant(glm::mat3(glm::vec3(P1.z, P1.x, 1), glm::vec3(P2.z, P2.x, 1), glm::vec3(P3.z, P3.x, 1)));
	double det3 = glm::determinant(glm::mat3(glm::vec3(P1.x, P1.y, 1), glm::vec3(P2.x, P2.y, 1), glm::vec3(P3.x, P3.y, 1)));
	double area = 0.5 * sqrt(det1 * det1 + det2 * det2  + det3 * det3);
	return area;
}

bool Mesh::lineIntersect(glm::vec3 const& p_start, glm::vec3 const& p_end, glm::vec3& normal, glm::vec3& intersect)
{
	//glm::vec3 V = glm::normalize(V0);
	glm::vec3 dir = p_end - p_start;
	if(glm::length(dir) < EPSILON)
		return false;
	glm::vec3 pos = p_start;
	float t = -1.0f;
	unsigned int n = face_normals.size();
	for(unsigned int i = 0; i < n; ++i)
	{
		glm::vec3 P1 = vertices[indices[i * 3 + 0]];
		glm::vec3 P2 = vertices[indices[i * 3 + 1]];
		glm::vec3 P3 = vertices[indices[i * 3 + 2]];

		if(P1 == P2 || P1 == P3 || P2 == P3)
			continue;

		glm::vec3 face_normal = -face_normals[i];
		//glm::vec3 face_normal = glm::normalize(glm::cross(P2 - P1, P3 - P1));
		if(glm::dot(face_normal, dir) < 0)
			face_normal = -face_normal;
		float denominator = glm::dot(face_normal, dir);
		if((denominator < EPSILON) && (denominator > -EPSILON))
			continue;
		else
		{
			t = glm::dot(face_normal, P1 - pos) / denominator;
			if(t < EPSILON || t > 1)
				continue;
			else
			{
				glm::vec3 P = pos + t * dir;
				float s = triangleArea(P1, P2, P3);
				float s1 = triangleArea(P, P2, P3) / s;
				float s2 = triangleArea(P, P3, P1) / s;
				float s3 = triangleArea(P, P1, P2) / s;
				// better solution?
				if((s1 > EPSILON) && (s1 < 1) && (s2 > EPSILON) && (s2 < 1) && (s3 > EPSILON) && (s3 < 1) && (s1 + s2 + s3 - 1 < EPSILON) && (s1 + s2 + s3 - 1 > -EPSILON))
				{
					normal = face_normal;
					intersect = P;
					return true;
				}
			}
		}
	}
	return false;
}