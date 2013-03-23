#include "particleSystem.h"
/////////////////////Particle///////////////////////////////////
particle::particle(){
		
		vel = glm::vec3(0.0);
		accel = glm::vec3(0.0);
		rest_density = 1.0;
		actual_density = 1.0;
		viscosity = 1;
		gas_constant = 1;
		temperature = 100;		
		color_interface = glm::vec3(1,1,1);
		color_surface =  glm::vec3(1,1,1);
	};

particle::particle(glm::vec3 position){

		pos = position;
		vel = glm::vec3(0.0);
		accel = glm::vec3(0.0);
		rest_density = 1.0;
		actual_density = 1.0;
		viscosity = 1;
		gas_constant = 1;
		temperature = 100;
		color_interface = glm::vec3(1,1,1);
		color_surface =  glm::vec3(1,1,1);
		
	};

void particle::Draw(const VBO& vbos, float r, int nSlice, int nStack){
		float phi   = 2.f * M_PI /(float)(nSlice-1);
		float theta = M_PI /(float)(nStack-1);

		int count = nSlice * nStack;

		std::vector<glm::vec3> vertices;
		std::vector<glm::vec3> normals;
		std::vector<glm::vec3> indices;
		std::vector<glm::vec3> colors;

		vertices.resize(count);
		colors.resize(count);
		normals.resize(count);
		indices.resize(count * 2);
		
		for(int i = 0; i < nStack; i++)
			for(int j = 0; j < nSlice; j++){

				float x = sin(j*phi) * cos(i*theta);
				float y = sin(j*phi) * sin(i*theta);
				float z = cos(j*phi);

				int index = i * nSlice + j;

				vertices[index] = r * glm::vec3(x, y, z) + pos;
				normals[ index] = glm::vec3(x, y, z);
				colors[  index] = glm::vec3(1,0,0);
		}

		for(int i = 0; i < nStack-1; i++)
			for(int j = 0; j < nSlice-1; j++){
				int index = (i*(nSlice-1)+j)*2;
				indices[index].x = i * nSlice + j;
				indices[index].y = i * nSlice + j + 1;
				indices[index].z = (i+1) * nSlice + j + 1;

				indices[index + 1].x = (i+1) * nSlice + j;
				indices[index + 1].y = i * nSlice + j;
				indices[index + 1].z = (i+1) * nSlice + j + 1;	
			}

	count *= 3;
	// position
    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_vbo);
    glBufferData(GL_ARRAY_BUFFER, count * sizeof(float), &vertices[0], GL_DYNAMIC_DRAW);
    // color
    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_cbo);
    glBufferData(GL_ARRAY_BUFFER, count * sizeof(float), &colors[0], GL_STATIC_DRAW);
    // normal
    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_nbo);
    glBufferData(GL_ARRAY_BUFFER, count * sizeof(float), &normals[0], GL_DYNAMIC_DRAW);

    // indices
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos.m_ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 2 * count * sizeof(unsigned int), &indices[0], GL_STATIC_DRAW);

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
    glDrawElements(GL_TRIANGLES, 2 * count, GL_UNSIGNED_INT, 0);

    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
};
//////////////////////ParticleSystem/////////////////////////////////
/////////////////////////////////////////////////////////////////////
particleSystem::particleSystem(){
}

particleSystem::particleSystem(int number){
	radius = 5.0;
	particles.resize(number*number*number);
	initParticles(number);

}

void particleSystem::initParticles(int number){
	for(int x = 0; x < number; x++)
		for(int y = 0; y < number; y++)
			for(int z = 0; z < number; z++){
				particle p(glm::vec3(x, y, z) * radius);
				particles.push_back(p);
			}
}

void particleSystem::Draw(const VBO& vbos){
	for (std::vector<particle>::iterator it = particles.begin() ; it != particles.end(); ++it)
		it->Draw(vbos, radius, 20, 20);
}

void particleSystem::LeapfrogIntegrate(){
	for (std::vector<particle>::iterator it = particles.begin() ; it != particles.end(); ++it){
	
	}
}