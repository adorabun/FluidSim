#include "particleSystem.h"
/////////////////////Particle///////////////////////////////////
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

		force = glm::vec3(0,9.8f,0);

		pos = position;
		vel = glm::vec3(0.0);

		rest_density = 1.0;
		actual_density = 1.0;
		viscosity = 1;
		gas_constant = 1;
		temperature = 100;
		color_interface = glm::vec3(1,0,0);
		color_surface =  glm::vec3(1,1,1);
		
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

	m_positions = p.m_positions;
	m_normals = p.m_normals;
	m_colors = p.m_colors;
	m_indices = p.m_indices;


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

	m_positions = p.m_positions;
	m_normals = p.m_normals;
	m_colors = p.m_colors;
	m_indices = p.m_indices;
    return *this;
}

void particle::Draw2(const VBO& vbos, float r, int nSlice, int nStack){
	m_positions.clear();
	m_normals.clear();
	m_colors.clear();
	m_indices.clear();

    glm::vec3 mat_color= color_interface ;

    glm::vec3 tnormal(0.0f, 1.0f, 0.0f), tpos;
	tpos = pos + r * tnormal;

    m_positions.push_back(tpos);
    m_normals.push_back(tnormal);
    m_colors.push_back(mat_color);

	float theta_z, theta_y, sin_z;
    float delta_y = 360.0f / nSlice, delta_z = 180.0f / nStack;
	//loop over the sphere
	for(theta_z = delta_z; theta_z < 179.99f; theta_z += delta_z)
	{
		for(theta_y = 0.0f; theta_y < 359.99f; theta_y += delta_y)
		{
			sin_z = sin(glm::radians(theta_z));
			
            tnormal.x = sin_z * cos(glm::radians(theta_y));
			tnormal.y = cos(glm::radians(theta_z));
			tnormal.z = -sin_z * sin(glm::radians(theta_y));

			tpos = pos + r * tnormal;

            m_positions.push_back(tpos);
            m_normals.push_back(tnormal);
            m_colors.push_back(mat_color);
		}
	}
	tnormal = glm::vec3(0.0f, -1.0f, 0.0f);
    tpos = pos + r * tnormal;

    m_positions.push_back(tpos);
    m_normals.push_back(tnormal);
    m_colors.push_back(mat_color);

	//indices
	unsigned int j = 0, k = 0;
	for(j = 0; j < nSlice - 1; ++j)
	{
		m_indices.push_back(0);
		m_indices.push_back(j + 1);
		m_indices.push_back(j + 2);
	}
	m_indices.push_back(0);
	m_indices.push_back(nSlice);
	m_indices.push_back(1);

	for(j = 0; j < nStack- 2; ++j)
	{
		for(k = 1 + nSlice * j; k < nSlice * (j + 1); ++k)
		{
			m_indices.push_back(k);
			m_indices.push_back(k + nSlice);
			m_indices.push_back(k + nSlice + 1);

			m_indices.push_back(k);
			m_indices.push_back(k + nSlice + 1);
			m_indices.push_back(k + 1);
		}
		m_indices.push_back(k);
		m_indices.push_back(k + nSlice);
		m_indices.push_back(k + 1);

		m_indices.push_back(k);
		m_indices.push_back(k + 1);
		m_indices.push_back(k + 1 - nSlice);
	}

    unsigned int bottom_id = (nStack - 1) * nSlice + 1;
    unsigned int offset = bottom_id - nSlice;
	for(j = 0; j < nSlice - 1; ++j)
	{
		m_indices.push_back(j + offset);
		m_indices.push_back(bottom_id);
		m_indices.push_back(j + offset + 1);
	}
	m_indices.push_back(bottom_id - 1);
	m_indices.push_back(bottom_id);
	m_indices.push_back(offset);

	if(m_indices.size() != 6 * (nStack - 1) * nSlice)
		printf("indices number not correct!\n");

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
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
}

void particle::Draw(const VBO& vbos, float r, int nSlice, int nStack){
		float phi   = 2.f * M_PI /(float)(nSlice-1);
		float theta = M_PI /(float)(nStack-1);

		m_positions.clear();
		m_normals.clear();
		m_colors.clear();
		m_indices.clear();
		
		for(int i = 0; i < nStack; i++)
			for(int j = 0; j < nSlice; j++){

				float x = sin( j*phi ) * cos( i*theta );
				float y = sin( j*phi ) * sin( i*theta );
				float z = cos( j*phi );

				m_positions.push_back( r * glm::vec3(x, y, z) + pos );
				 m_normals.push_back(     glm::vec3(x, y, z)       );
				  m_colors.push_back(     glm::vec3(1, 0, 0)       );
		}

		for(int i = 0; i < nStack-1; i++)
			for(int j = 0; j < nSlice-1; j++){
				m_indices.push_back(     i * nSlice + j     );
				m_indices.push_back(     i * nSlice + j + 1 );
				m_indices.push_back( (i+1) * nSlice + j + 1 );

				m_indices.push_back( (i+1) * nSlice + j    );
				m_indices.push_back(     i * nSlice + j    );
				m_indices.push_back( (i+1) * nSlice + j + 1);
			}

	
	// position
    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_vbo);
	glBufferData(GL_ARRAY_BUFFER, 3 * m_positions.size() * sizeof(float), &m_positions[0], GL_DYNAMIC_DRAW);
    // color
    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_cbo);
    glBufferData(GL_ARRAY_BUFFER, 3 * m_colors.size() * sizeof(float), &m_colors[0], GL_STATIC_DRAW);
    // normal
    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_nbo);
    glBufferData(GL_ARRAY_BUFFER, 3 * m_normals.size() * sizeof(float), &m_normals[0], GL_DYNAMIC_DRAW);

    // indices
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos.m_ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_indices.size() * sizeof(unsigned short), &m_indices[0], GL_STATIC_DRAW);

	//activate our three kinds of information
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
	//draw the elements
    glDrawElements(GL_TRIANGLES, m_indices.size(), GL_UNSIGNED_INT, 0);

	//shut off the information since we're done drawing
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);

   // glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}


//////////////////////ParticleSystem/////////////////////////////////
/////////////////////////////////////////////////////////////////////
particleSystem::particleSystem(){
}

particleSystem::particleSystem(int number){
	radius = 0.15;
	initParticles(number);

}

void particleSystem::initParticles(int number){
	for(int x = 0; x < number; x++)
		for(int y = 0; y < number; y++)
			for(int z = 0; z < number; z++){
				particle p(glm::vec3(x, y, z) * radius * 2.f);
				particles.push_back(p);
			}
}

void particleSystem::Draw(const VBO& vbos){
	LeapfrogIntegrate(0.01);
	for (std::vector<particle>::iterator it = particles.begin() ; it != particles.end(); ++it)
		it->Draw2(vbos, radius, 6, 6);
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