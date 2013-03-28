#include "particleSystem.h"
/////////////////////Particle///////////////////////////////////
#define nSlice 6
#define nStack 6
#define radius 0.15f

particle::particle(){
		
		mass = 1.f;
		vel = glm::vec3(0.0);
		
		rest_density = 1.0;
		actual_density = 1.0;
		viscosity = 1;
		gas_constant = 1;
		temperature = 100;		
		color_interface = glm::vec3(1,1,1);
		color_surface =  glm::vec3(1,1,1);
	}

particle::particle(glm::vec3 position){
		
		mass = 1.f;

		force = glm::vec3(0,-9.8f,0);

		pos = position;
		vel = glm::vec3(0.0);

		rest_density = 1.0;
		actual_density = 1.0;
		viscosity = 1;
		gas_constant = 1;
		temperature = 100;
		color_interface = glm::vec3(0,1,0);
		color_surface =  glm::vec3(1,0,1);

		
	}

particle::particle(const particle& p) 
{
	mass = p.mass;
	pos = p.pos;
	vel = p.vel;
	force = p.force;
	rest_density = p.rest_density;
	actual_density = p.actual_density;

	viscosity = p.viscosity;
	gas_constant = p.gas_constant;
	temperature = p.temperature;

	color_interface = p.color_interface;
	color_surface = p.color_surface;

}

particle& particle::operator=(const particle& p)
{
    if (&p == this) return *this;

	mass = p.mass;
	pos = p.pos;
	vel = p.vel;
	force = p.force;
	rest_density = p.rest_density;
	actual_density = p.actual_density;

	viscosity = p.viscosity;
	gas_constant = p.gas_constant;
	temperature = p.temperature;

	color_interface = p.color_interface;
	color_surface = p.color_surface;

    return *this;
}



//////////////////////ParticleSystem/////////////////////////////////
/////////////////////////////////////////////////////////////////////
particleSystem::particleSystem(){
}

particleSystem::particleSystem(int number){
	
	initParticles(number);

}

void particleSystem::initParticles(int number){

	float stepsize = 2.f * radius;

	particles.resize(number * number * number);

	for(int x = 0; x < number; x++)
		for(int y = 0; y < number; y++)
			for(int z = 0; z < number; z++){

				particle p(glm::vec3(x, y, z) * stepsize);

				particles[ x*number*number + y*number + z] = p;
			}

	initSphere();
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

void particleSystem::Draw(const VBO& vbos){
	LeapfrogIntegrate(0.01);
	int index;

	for (std::vector<particle>::iterator it = particles.begin() ; it != particles.end(); ++it){

		for(int i = 0; i < nStack; i++)
			for(int j = 0; j < nSlice; j++){
				index = i*nSlice+j;
				m_positions[index] += (it->pos);
				m_colors[index] = it->color_surface;
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

void particleSystem::LeapfrogIntegrate(float dt){
	float halfdt = 0.5 * dt;
	particleGrid target = particles;// target is a copy!
	for (int i=0; i < target.size(); i++){
		target[i].pos = particles[i].pos + particles[i].vel * dt + halfdt * dt * particles[i].force * 1.f / particles[i].mass;
	}

	for (int i=0; i < target.size(); i++){
		particles[i].vel += halfdt * (target[i].force  + particles[i].force) * 1.f / particles[i].mass;
		particles[i].pos = target[i].pos;
	}

}