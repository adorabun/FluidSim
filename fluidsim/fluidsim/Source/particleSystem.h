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

	void Draw(const VBO& vbos, float r, int nSlice, int nStack);

	float mass;
	glm::vec3 pos;
	glm::vec3 vel;
	glm::vec3 accel;
	 
	float rest_density;
	float actual_density;

	float viscosity;

	float gas_constant;

	glm::vec3 color_interface;
	glm::vec3 color_surface;

	float temperature;

protected:
};

class particleSystem
{
public:
	particleSystem();
	particleSystem(int number);
	void Draw(const VBO& vbos);
	void LeapfrogIntegrate();

private:
	void initParticles(int number);
	std::vector<particle> particles;
	float radius;

	
	
};

#endif

