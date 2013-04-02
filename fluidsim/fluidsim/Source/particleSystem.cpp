#include "particleSystem.h"
#include <cmath>
#include <fstream>
#include <iostream>

#define EPSILON 0.000001f
#define SMOOTH_CORE_RADIUS 1.f

float particleSystem::nSlice = 6;
float particleSystem::nStack = 6;
float particleSystem::radius = 0.15;

/////////////////////Particle///////////////////////////////////
particle::particle(){

	}

particle::particle(glm::vec3 position){
		
		mass = 19.683f;
		force = glm::vec3(0.f);

		pos = position;
		vel = glm::vec3(0.0);

		rest_density = 1000.f;
		actual_density = rest_density;

		viscosity_coef = 10.f;
		gas_constant = 3.f;
		
		temperature = 500;

		color_interface = 1.f;
		color_surface = 1.f;
		
	


	}


//////////////////////ParticleSystem/////////////////////////////////
particleSystem::particleSystem(){

}

particleSystem::particleSystem(int number){
	xstart = 0.f;
	ystart = 0.f;
	zstart = 0.f;
	xend = 4.f;
	yend = 4.f;
	zend = 4.f;

	initParticles(number);
	initSphere();
	//gridcells.resize(4,4,4, particles);
}

void particleSystem::initParticles(int number){

	float stepsize = 2.f * radius;

	
	particles.resize(number * number * number);
	
	for(int x = 0; x < number; x++)
		for(int y = 0; y < number; y++)
			for(int z = 0; z < number; z++){

				particle p(glm::vec3(x, y, z) * stepsize + glm::vec3(radius));

				particles[x*number*number + y*number + z] = p;
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
/*neighbor search
each cells has index of particles
*/
//do neighbor search
//compute density 
//compute force
//integrate
//scene interaction
//visualization
void particleSystem::LeapfrogIntegrate(float dt){
	float halfdt = 0.5f * dt;
	particleGrid target = particles;// target is a copy!
	particleGrid& source = particles;//source is a ptr!

	for (int i=0; i < target.size(); i++){

		//what should i do if actual_density is 0?
		
		/*if(source[i].actual_density == 0.f)
			continue;*/
		target[i].pos = source[i].pos + source[i].vel * dt 
						+ halfdt * dt * source[i].force / source[i].actual_density;
		
	}

	//particleGrid nghrs;
	glm::vec3 collision_normal = glm::vec3(0.f);
	//calculate actual density 
	for (int i=0; i < target.size(); i++){
		//if( !checkIfOutOfBoundry(target[i]) ){
		//	nghrs = gridcells.getNeighbors(target[i]);
		//	target[i].actual_density = computeDensity(nghrs, target[i]);
		//	
		//	target[i].pressure = target[i].gas_constant * (target[i].actual_density - target[i].rest_density);
		//}else{
			target[i].actual_density = computeDensity(target, target[i]);
		
			target[i].pressure = target[i].gas_constant * (target[i].actual_density - target[i].rest_density);
		//}
	}

	
	for (int i=0; i < target.size(); i++){
			//if( !checkIfOutOfBoundry(target[i]) ){
			//	nghrs = gridcells.getNeighbors(target[i]);
			//	target[i].force = computeForce(nghrs, target[i]);
			//}else

		//what should i do if actual_density is 0?
		/*if (target[i].actual_density == 0.f)
			continue;
		*/target[i].force = computeForce(target, target[i]);
		
	}

	
	for (int i=0; i < target.size(); i++){
		if( CollisionDectection2(target[i], collision_normal) ){
			//glm::vec3 vn = (source[i].vel * collision_normal) * collision_normal;//decompose v along normal
			//glm::vec3 vt = source[i].vel - vt;
			//source[i].vel = 0.9f * vt - 0.5f * vn;//flip normal direction speed
			//
			source[i].vel.y *= -1.f;
		}else{
			//what should i do if actual_density is 0?
		/*	if(source[i].actual_density == 0.f && target[i].actual_density == 0.f)
				continue;
			else if(source[i].actual_density == 0.f)
				source[i].vel += halfdt * target[i].force/target[i].actual_density;
			else if (target[i].actual_density == 0.f)
				source[i].vel += halfdt * source[i].force /source[i].actual_density;
			else*/
			source[i].vel += halfdt * (target[i].force/target[i].actual_density  + source[i].force /source[i].actual_density);
			source[i].pos = target[i].pos;
			source[i].force = target[i].force;
			source[i].actual_density = target[i].actual_density;
		}
	}

	//gridcells.refillGrid(source);

}

bool particleSystem::checkIfOutOfBoundry(particle p){
	glm::vec3 pos = p.pos;
	if(pos.x < xstart || pos.x > xend || pos.y < ystart || pos.y > yend || pos.z < zstart || pos.z > zend)
		return true;
	return false;
}

bool particleSystem::CollisionDectection(particle p, glm::vec3& n){
	n = glm::vec3(0.f);
	glm::vec3 pos = p.pos;

	if(pos.x < xstart + EPSILON)
		n.x = 1.f;
	else if(pos.x > xend - EPSILON)
		n.x = -1.f;
	
	if(pos.y < ystart + EPSILON)
		n.y = 1.f;
	else if(pos.y > yend - EPSILON)
		n.y = -1.f;

	if(pos.z < zstart + EPSILON)
		n.z = 1.f;
	else if(pos.z > zend - EPSILON)
		n.z = -1.f;

	if(n == glm::vec3(0.f))
		return false;

	glm::normalize(n);
	return true;
}

bool particleSystem::CollisionDectection2(particle p, glm::vec3& n){
	n = glm::vec3(0.f);
	glm::vec3 pos = p.pos;

	
	if(pos.y <= ystart)
		n.y = 1.f;
	else if(pos.y >= yend)
		n.y = -1.f;

	if(n == glm::vec3(0.f))
		return false;

	glm::normalize(n);
	return true;
}
///////////////////////////////smoothing kernels//////////////////////////////////////
inline float poly6Kernel(glm::vec3 r, float h){
	float rLen = glm::length(r);

	if(rLen <= h){
		float t = 315.f * pow(h*h-rLen*rLen, 3) / (64.f * M_PI * pow(h,9) );
		return t;
	}
	return 0.f;
}

inline glm::vec3 poly6KernelGradient(glm::vec3 r, float h){
	float rLen = glm::length(r);

	if(rLen <= h){
		float t =  945.f  * pow(h*h - rLen*rLen, 2) / (32.f * M_PI * pow(h,9) );
		return t*r;
	}
		
	return glm::vec3(0.f);

}

inline float poly6KernelLaplacian(glm::vec3 r, float h){
	float rLen = glm::length(r);

	if(rLen <= h)
		return 945.f * (h*h - rLen*rLen) * (7.f*rLen*rLen - 3*h*h) / (32.f * M_PI * pow(h,9) );
		
	return 0.f;
}

inline float spikyKernel(glm::vec3 r, float h){
	float rLen = glm::length(r);

	if(rLen <= h)
		return 15.f * pow(h - rLen, 3) / ( M_PI * pow(h,6) ) ;
	return 0.f;
}

inline glm::vec3 spikyKernelGradient(glm::vec3 r, float h){

	float rLen = glm::length(r);

	if(rLen > 0 && rLen <= h){
		float t = - 45.f * ( ( h*h + rLen*rLen ) / rLen - 2.f * h ) / (M_PI * pow(h, 6));
		return t * r;
	}
	return glm::vec3(0.f);
}

inline float viscosityKernel(glm::vec3 r, float h){
	float rLen = glm::length(r);

	if(rLen <= h)
		return 15.f * ( - rLen*rLen*rLen/(2.f*h*h*h) + rLen*rLen/(h*h) + h/(2.f*rLen) - 1 ) / ( 2 * M_PI * pow( h, 3 ) ) ;
	return 0.f;
}

inline float viscosityKernelLaplacian(glm::vec3 r, float h){
	float rLen = glm::length(r);
	if(rLen <= h)
		return 45.f * ( h - rLen ) / ( M_PI + pow(h, 6) ) ;
	return 0.f;
}

///////////////////////////////computation//////////////////////////////////////
glm::vec3 particleSystem::computeForce(const particleGrid& ps, particle pi){
	glm::vec3 f_pressure(0.f);
	glm::vec3 f_viscosity(0.f);
	glm::vec3 f_surfaceTension(0.f);
	glm::vec3 f_gravity = glm::vec3(0,-9.8f,0) * pi.actual_density;
	
	float massOverDensity;
	glm::vec3 r;

	float tension_coeff = 1.f;
	glm::vec3 Cs_normal = glm::vec3(0.f);
	float Cs_Laplacian = 0.f;

	for (int j=0; j< ps.size(); j++){
		massOverDensity = ps[j].mass / ps[j].actual_density;

		r = pi.pos - ps[j].pos;

		
		f_pressure -=  massOverDensity * ( pi.pressure + ps[j].pressure ) * 0.5f * spikyKernelGradient(r, SMOOTH_CORE_RADIUS);
		f_viscosity  += massOverDensity * ( ps[j].vel - pi.vel ) * viscosityKernelLaplacian(r, SMOOTH_CORE_RADIUS);

		Cs_normal += massOverDensity * poly6KernelGradient(r, SMOOTH_CORE_RADIUS);
		Cs_Laplacian += massOverDensity * poly6KernelLaplacian(r, SMOOTH_CORE_RADIUS);
		
	}

	f_viscosity *= pi.viscosity_coef;

	
	float sCs_normal_len = glm::length(Cs_normal);
	if(sCs_normal_len > 0.8f){
		float curvature =  - Cs_Laplacian / sCs_normal_len;
		f_surfaceTension = tension_coeff * curvature * Cs_normal;
	}

	return f_pressure + f_viscosity + f_surfaceTension + f_gravity;
	//return f_gravity;
}

float particleSystem::computeDensity(const particleGrid& ps, particle pi){

	//resrt
	float rho = 0.f;

	//accumulate
	for (int j=0; j< ps.size(); j++)
	{
		rho  += ps[j].mass * poly6Kernel(pi.pos - ps[j].pos, SMOOTH_CORE_RADIUS);
	}

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

void particleSystem::outputCenter(int& i_frame, char* s_file)
{
	std::ofstream io_out;
	if(i_frame == 0)
		io_out.open(s_file, std::ios::ate);
	else
		io_out.open(s_file, std::ios::app);
	io_out<<i_frame<<" ";
	glm::vec3 center = glm::vec3(0.0,0.0,0.0);
	int n = 0;
	std::cout<<particles.size();
	for (std::vector<particle>::iterator it = particles.begin() ; it != particles.end(); ++it)
	{
		center = it->pos;
		io_out<<center.x<<" "<<center.y<<" "<<center.z<<" ";
		n++;
		if(n > 3)
			break;
	}
	io_out<<std::endl;
	i_frame++;
}
///////////////////////////////Grid//////////////////////////////////////
//void particleSystem::Grid::resize(int x, int y, int z, const particleGrid& ps){
//	dim = glm::vec3(x,y,z);
//	GridData.resize(dim.x * dim.y * dim.z);
//	refillGrid(ps);
//}
//
//void particleSystem::Grid::refillGrid(const particleGrid& ps){
//	for(int i=0; i < GridData.size(); i++){
//		GridData[i].clear();
//	}
//
//	for(int i=0; i < ps.size(); i++){
//		pushParticle(ps[i]);
//	}
//}
//
//void particleSystem::Grid::pushParticle(const particle& p){
//	particle np = p;
//	int index = positionToIndex(np.pos);
//	GridData[index].push_back(np);
//}
//particleSystem::particleGrid particleSystem::Grid::getNeighbors(const particle& p){
//	int index = positionToIndex(p.pos);
//	return GridData[index];
//}
//
//int particleSystem::Grid::positionToIndex(glm::vec3 p){
//	glm::vec3 index(0);
//
//	index.x = (int) (p.x/SMOOTH_CORE_RADIUS);
//	index.y = (int) (p.y/SMOOTH_CORE_RADIUS);
//	index.z = (int) (p.z/SMOOTH_CORE_RADIUS);
//
//	return gridToVec(index);
//}
//int particleSystem::Grid::gridToVec(glm::vec3 index){
//
//	return dim.x * dim.y * index.z + dim.x * index.y + index.x;
//
//}

