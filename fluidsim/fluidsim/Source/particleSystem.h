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

	int id;
	float mass;
	
	glm::vec3 pos;
	glm::vec3 vel;
	glm::vec3 force;
	 
	float rest_density;
	float actual_density;

	float viscosity_coef;//mu
	float gas_constant;//k

	float temperature;

	float color_interface;//Ci
	float color_surface;//Cs

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
	static float nSlice;
	static float nStack;
	static float radius;
	static float tension_coeff;//sigma
	static float surfaceThreshold;//l
	static float xstart;
	static float ystart;
	static float zstart;
	static float xend;
	static float yend;
	static float zend;

	void outputCenter(int& i_frame, char* s_file);

private:
	void initParticles(int number);
	void initSphere();

	
	void computeForce(const particleGrid& ps, particle& pi);
	void computeDensity(const particleGrid& ps, particle& pi);
	bool checkIfOutOfBoundry(const particle& p);
	bool CollisionDectection(const particle& p, glm::vec3& n);

private:
	particleGrid particles;

	class Grid{
	public:
		void resize(int x, int y, int z, const particleGrid& ps);
		void refillGrid(const particleGrid& ps);

		void pushParticle(const particle& p);
		void getNeighbors(const particleGrid& ps, const particle& p, particleGrid& des);
		

		glm::vec3 positionToGridIndex(glm::vec3 p);
		int gridIndexToVecIndex(glm::vec3 index);
		int positionToVecIndex(glm::vec3 p);
		bool IfWithinBoundry(glm::vec3 gridIndex);
		
		glm::vec3 dim;

		std::vector<std::vector<int>> GridData;
	};

	Grid gridcells;

	std::vector<glm::vec3> m_positions;
	std::vector<glm::vec3> m_normals;
	std::vector<glm::vec3> m_colors;
	std::vector<unsigned short> m_indices;

	
};

#endif

