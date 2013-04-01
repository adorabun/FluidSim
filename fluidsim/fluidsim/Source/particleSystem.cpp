#include "particleSystem.h"
#include <cmath>

float particleSystem::nSlice = 6;
float particleSystem::nStack = 6;
float particleSystem::radius = 0.15;

/////////////////////Particle///////////////////////////////////
particle::particle(){

	}

particle::particle(glm::vec3 position){
		
		mass = 1.f;
		force = glm::vec3(0.f);

		pos = position;
		vel = glm::vec3(0.0);

		rest_density = 1.f;
		actual_density = 1.f;

		viscosity_coef = 1.f;
		gas_constant = 1.f;
		
		temperature = 100;

		color_interface = 1.f;
		color_surface = 1.f;


	}


//////////////////////ParticleSystem/////////////////////////////////
particleSystem::particleSystem(){

}

particleSystem::particleSystem(int number){
	xstart = -3.0;
	ystart = -3.0;
	zstart = -3.0;
	xend = 9.0;
	yend = 9.0;
	zend = 9.0;

	initParticles(number);
	initSphere();
}

void particleSystem::initParticles(int number){

	float stepsize = 2.f * radius;

	particles.resize(number * number * number);

	for(int x = 0; x < number; x++)
		for(int y = 0; y < number; y++)
			for(int z = 0; z < number; z++){

				particle p(glm::vec3(x, y, z) * stepsize + glm::vec3(radius));

				particles[ x*number*number + y*number + z] = p;
			}
	
	for (int i=0; i < particles.size(); i++){
		particles[i].actual_density = computeDensity(particles, i);
		particles[i].pressure = particles[i].gas_constant * (particles[i].actual_density - particles[i].rest_density);
	}
}

void particleSystem::initSphere(){

	float phi   = 2.f * M_PI /(float)(nSlice-1);
	float theta = M_PI /(float)(nStack-1);


	int count = nSlice*nStack;

	m_positions.resize( count );
	m_normals.resize( count );
	m_colors.resize( count );
	m_indices.resize( count * 6 );

	for(int i = 0; i < nStack; i++)
		for(int j = 0; j < nSlice; j++){

			float x = sin(j*phi) * cos(i*theta);
			float y = sin(j*phi) * sin(i*theta);
			float z = cos(j*phi);

			int index = i*nSlice+j;
			m_positions[index] = glm::vec3(x, y, z) * radius;
			m_normals[index] = glm::vec3(x, y, z);
		}

	for(int i = 0; i < nStack-1; i++)
		for(int j = 0; j < nSlice-1; j++){
			int index = (i*(nSlice-1)+j)*6;
			m_indices[index    ] = i * nSlice + j;
			m_indices[index + 1] = i * nSlice + j + 1;
			m_indices[index + 2] = (i+1) * nSlice + j + 1;

			m_indices[index + 3] = (i+1) * nSlice + j;
			m_indices[index + 4] = i * nSlice + j;
			m_indices[index + 5] = (i+1) * nSlice + j + 1;	
	}

}


//how to handle boundary case?
//how to speed up neighboring search?
//make a particle grid
void particleSystem::LeapfrogIntegrate(float dt){
	float halfdt = 0.5f * dt;
	particleGrid target = particles;// target is a copy!
	particleGrid& source = particles;//source is a ptr!

	for (int i=0; i < target.size(); i++){
		target[i].pos = source[i].pos + source[i].vel * dt 
						+ halfdt * dt * source[i].force / source[i].actual_density;
	}

	//calculate actual density 
	
	for (int i=0; i < target.size(); i++){
		target[i].actual_density = computeDensity(target, i);
		target[i].pressure = target[i].gas_constant * (target[i].actual_density - target[i].rest_density);
	}

	for (int i=0; i < target.size(); i++){
		target[i].force = computeForce(target, i);
	}

	for (int i=0; i < target.size(); i++){
		source[i].vel += halfdt * (target[i].force/target[i].actual_density  + source[i].force /source[i].actual_density);
		source[i].pos = target[i].pos;
	}

}

///////////////////////////////smoothing kernels//////////////////////////////////////
inline float poly6Kernel(glm::vec3 r, float h){
	float rLen = glm::length(r);

	if(rLen >= 0.f || rLen <= h)
		return 315.f * pow(h*h-rLen*rLen, 3) / (64.f * M_PI * pow(h,9) ) ;
	return 0.f;
}

inline glm::vec3 poly6KernelGradient(glm::vec3 r, float h){
	float rLen = glm::length(r);

	if(rLen >= 0.f || rLen <= h){
		float t =  945.f  * pow(h*h - rLen*rLen, 2) / (32.f * M_PI * pow(h,9) );
		return t*r;
	}
		
	return glm::vec3(0.f);

}

inline float poly6KernelLaplacian(glm::vec3 r, float h){
	float rLen = glm::length(r);

	if(rLen >= 0.f || rLen <= h)
		return 945.f * (h*h - rLen*rLen) * (7.f*rLen*rLen - 3*h*h) / (32.f * M_PI * pow(h,9) );
		
	return 0.f;
}

inline float spikyKernel(glm::vec3 r, float h){
	float rLen = glm::length(r);

	if(rLen >= 0.f || rLen <= h)
		return 15.f * pow(h - rLen, 3) / ( M_PI * pow(h,6) ) ;
	return 0.f;
}

inline glm::vec3 spikyKernelGradient(glm::vec3 r, float h){

	float rLen = glm::length(r);
	if(rLen >= 0.f || rLen <= h){
		float t = 45.f * ( ( h*h + rLen*rLen ) / rLen - 2.f * h ) / (M_PI * pow(h, 6));
		return t * r;
	}
	return glm::vec3(0.f);
}

inline float viscosityKernel(glm::vec3 r, float h){
	float rLen = glm::length(r);

	if(rLen >= 0.f || rLen <= h)
		return 15.f * ( - rLen*rLen*rLen/(2.f*h*h*h) + rLen*rLen/(h*h) + h/(2.f*rLen) - 1 ) / ( 2 * M_PI * pow( h, 3 ) ) ;
	return 0.f;
}

inline float viscosityKernelLaplacian(glm::vec3 r, float h){
	float rLen = glm::length(r);
	if(rLen >= 0.f || rLen <= h)
		return 45.f * ( h - rLen ) / ( M_PI + pow(h, 6) ) ;
	return 0.f;
}

///////////////////////////////computation//////////////////////////////////////
glm::vec3 particleSystem::computeForce(const particleGrid& ps, int index){
	glm::vec3 f_pressure(0.f);
	glm::vec3 f_viscosity(0.f);
	glm::vec3 f_surfaceTension(0.f);
	glm::vec3 f_gravity = glm::vec3(0,-9.8f,0);
	
	float massTimesInvDensity;
	glm::vec3 r;

	float tension_coeff = 1.f;
	glm::vec3 Cs_normal = glm::vec3(0.f);
	float Cs_Laplacian = 0.f;

	for (int i=0; i< ps.size(); i++){
		massTimesInvDensity = ps[i].mass / ps[i].actual_density;
		r = ps[index].pos - ps[i].pos;

		
		f_pressure -=  massTimesInvDensity * ( ps[index].pressure + ps[i].pressure ) * 0.5f * spikyKernelGradient(r, 1.f);
		f_viscosity  += massTimesInvDensity* ( ps[index].vel - ps[i].vel ) * viscosityKernelLaplacian(r, 1.f);

		Cs_normal += massTimesInvDensity * poly6KernelGradient(r, 1.f);
		Cs_Laplacian += massTimesInvDensity * poly6KernelLaplacian(r, 1.f);
		
	}

	f_viscosity *= ps[index].viscosity_coef;
	f_surfaceTension = - tension_coeff * Cs_Laplacian * glm::normalize(Cs_normal);

	//return f_pressure + f_viscosity + f_surfaceTension + f_gravity;
	return f_gravity;
}

float particleSystem::computeDensity(const particleGrid& ps, int index){

	//resrt
	float rho = 0.f;

	//accumulate
	for (int i=0; i< ps.size(); i++)
		rho  += ps[i].mass * poly6Kernel(ps[index].pos - ps[i].pos, 1.f);
	
	return rho;
		
}

///////////////////////////////draw related//////////////////////////////////////
void particleSystem::Draw(const VBO& vbos){
	
	LeapfrogIntegrate(0.01f);
	
	int index;

	for (std::vector<particle>::iterator it = particles.begin() ; it != particles.end(); ++it){

		for(int i = 0; i < nStack; i++)
			for(int j = 0; j < nSlice; j++){
				index = i*nSlice+j;
				m_positions[index] += (it->pos);
				m_colors[index] = glm::vec3(0.f,0.f, 1.f);
			}
		 // position
		glBindBuffer(GL_ARRAY_BUFFER, vbos.m_vbo);
		glBufferData(GL_ARRAY_BUFFER, 3 * m_positions.size() * sizeof(float), &m_positions[0], GL_STREAM_DRAW);

		// color
		glBindBuffer(GL_ARRAY_BUFFER, vbos.m_cbo);
		glBufferData(GL_ARRAY_BUFFER, 3 * m_colors.size() * sizeof(float), &m_colors[0], GL_STREAM_DRAW);

		// normal
		glBindBuffer(GL_ARRAY_BUFFER, vbos.m_nbo);
		glBufferData(GL_ARRAY_BUFFER, 3 * m_normals.size() * sizeof(float), &m_normals[0], GL_STREAM_DRAW);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos.m_ibo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_indices.size() * sizeof(unsigned short), &m_indices[0], GL_STATIC_DRAW);

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
		glDrawElements(GL_TRIANGLES, m_indices.size(), GL_UNSIGNED_SHORT, 0);//GL_UNSIGNED_INT

		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(1);
		glDisableVertexAttribArray(2);

		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

		for(int i = 0; i < nStack; i++)
			for(int j = 0; j < nSlice; j++){
				m_positions[i*nSlice+j] -= (it->pos);
		}

	}
}



void particleSystem::drawWireGrid()
{
   // Display grid in light grey, draw top & bottom

	  glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
      glDisable(GL_LIGHTING);
      glColor3f(0.25, 0.25, 0.25);

      glBegin(GL_LINES);
	  
	  glVertex3d(xstart, ystart, zstart);
	  glVertex3d(xend, ystart, zstart);

	  glVertex3d(xstart, yend, zstart);
	  glVertex3d(xend, yend, zstart);

	  glVertex3d(xstart, ystart, zend);
	  glVertex3d(xend, ystart, zend);

	  glVertex3d(xstart, yend, zend);
	  glVertex3d(xend, yend, zend);

	  glVertex3d(xstart, ystart, zstart);
	  glVertex3d(xstart, ystart, zend);

	   glVertex3d(xend, ystart, zstart);
	   glVertex3d(xend, ystart, zend);


	   glVertex3d(xstart, yend, zstart);
	   glVertex3d(xstart, yend, zend);

	  glVertex3d(xend, yend, zend);
	  glVertex3d(xend, yend, zstart);

      glVertex3d(xstart, ystart, zstart);
      glVertex3d(xstart, yend, zstart);

      glVertex3d(xend, ystart, zstart);
      glVertex3d(xend, yend, zstart);

      glVertex3d(xstart, ystart, zend);
      glVertex3d(xstart, yend, zend);

      glVertex3d(xend, ystart, zend);
      glVertex3d(xend, yend, zend);
      glEnd();
   glPopAttrib();

   glEnd();
}