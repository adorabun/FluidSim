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

	float mass;
	
	glm::vec3 pos;
	glm::vec3 vel;
	glm::vec3 force;
	 
	float rest_density;
	float actual_density;

	float viscosity_coef;

	float gas_constant;

	

	float temperature;

	glm::vec3 color_interface;
	glm::vec3 color_surface;

	float pressure;
};

class particleSystem
{
public:
	particleSystem();
	particleSystem(int number);

	void Draw(const VBO& vbos);
	void LeapfrogIntegrate(float dt);
	void drawWireGrid();

	typedef std::vector<particle> particleGrid;

private:
	void initParticles(int number);
	void initSphere();

	
	glm::vec3 computeForce(const particleGrid& ps, int index);

	//glm::vec3 computeColor(const particleGrid& ps, int index);
	float computeDensity(const particleGrid& ps, int index);
	
private:
	particleGrid particles;

	std::vector<glm::vec3> m_positions;
	std::vector<glm::vec3> m_normals;
	std::vector<glm::vec3> m_colors;
	std::vector<unsigned short> m_indices;
	
	double xstart;
	double ystart;
	double zstart;
	double xend;
	double yend;
	double zend;
};

#endif

