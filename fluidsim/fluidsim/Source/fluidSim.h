#ifndef FLUID_SIM_H_
#define FLUID_SIM_H_

#include "particleSystem.h"

class fluidSim{
public:
	fluidSim();
	~fluidSim();
	void SPH_solver();
	void update();
private:
	//simulation
	void initParticles();
	void LeapfrogIntegrate();
	void computePressure();
	void computeViscosity();
	void computeSurfaceTension();
	void computeExternalForces();
	
private:
	particleSystem ps;
};

#endif








