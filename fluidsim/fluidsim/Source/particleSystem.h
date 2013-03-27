#ifndef PARTICLE_SYSTEM_H_INCLUDED
#define PARTICLE_SYSTEM_H_INCLUDED

#include <vector>

#include "openGL_headers.h"
#include "math_headers.h"

#define M_PI 3.1415926

class particle{
public:
	particle();

	particle(glm::vec3 position);

	particle(const particle& p);
    particle& operator=(const particle& p);

	void Draw(const VBO& vbos, float r, int nSlice, int nStack);
	void Draw2(const VBO& vbos, float r, int nSlice, int nStack);
	
	float mass;
	
	glm::vec3 pos;
	glm::vec3 vel;

	glm::vec3 force;
	 
	float rest_density;
	float actual_density;

	float viscosity;

	float gas_constant;

	glm::vec3 color_interface;
	glm::vec3 color_surface;

	float temperature;

protected:
		std::vector<glm::vec3> m_positions;
		std::vector<glm::vec3> m_normals;
		std::vector<glm::vec3> m_colors;
		std::vector<unsigned short> m_indices;
};

class particleSystem
{
public:
	particleSystem();
	particleSystem(int number);

	void Draw(const VBO& vbos);
	void LeapfrogIntegrate(float dt);

	typedef std::vector<particle> particleGrid;

private:
	void initParticles(int number);
	
	particleGrid particles;
	float radius;

	
	
};

#endif

