#ifndef PARTICLE_SYSTEM_H_INCLUDED
#define PARTICLE_SYSTEM_H_INCLUDED

#include <vector>
#include <map>


#include "openGL_headers.h"
#include "math_headers.h"

#define M_PI 3.1415926

class particle{
public:
	particle();

	particle(glm::vec3 position, float rho);

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

	std::vector<particle*> ngbrs;
	
};

struct Cell{
	int frameID;
	std::vector<particle*> ps;
	Cell():frameID(-1){};
};

class SpaceGrid{
	public:
		void pushParticle(particle& pt, int frameID);
		void getNeighbors(particle& pt, int frameID);
		glm::vec3 dim;

	protected:
		
		glm::vec3 positionToGridIndex(glm::vec3 p);
		int gridIndexToVecIndex(glm::vec3 index);
		int positionToVecIndex(glm::vec3 p);
		bool IfWithinBoundry(glm::vec3 gridIndex);
		
		std::map<int, Cell> mymap;

};

class particleSystem
{
public:
	particleSystem();
	particleSystem(int numberX, int numberY, int numberZ);

	int frameCount;
	void Draw(const VBO& vbos);
	void LeapfrogIntegrate(float dt);
	void drawWireGrid();

	static int nSlice;
	static int nStack;
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
	void initParticles(int numberX, int numberY, int numberZ);
	void initSphere();
	void initCube();
	
	void GenerateParticles(int numberX, int numberY, int numberZ);
	void computeForce(particle& pi);
	void computeDensity(particle& pi);
	bool checkIfOutOfBoundry(const particle& p);
	bool CollisionDectection(particle& p, glm::vec3& n);

private:
	std::vector<particle> particles;


	SpaceGrid mygrid;

	std::vector<glm::vec3> m_positions;
	std::vector<glm::vec3> m_normals;
	std::vector<glm::vec3> m_colors;
	std::vector<unsigned short> m_indices;

	
};

#endif

