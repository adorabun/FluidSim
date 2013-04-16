#include "particleSystem.h"
#include <cmath>
#include <fstream>
#include <iostream>

#define EPSILON 0.000001f
#define SMOOTH_CORE_RADIUS 1.f //three times the average distance


int particleSystem::nSlice = 6;
int particleSystem::nStack = 6;
float particleSystem::radius = 0.15f;
float particleSystem::tension_coeff = 50.f;//sigma
float particleSystem::surfaceThreshold = 0.5f;//l need to figure out value by printing

float particleSystem::xstart = 0.f;

float particleSystem::ystart = 0.f;
float particleSystem::zstart = 0.f;
float particleSystem::xend = 4.0f;
float particleSystem::yend = 10.0f;
float particleSystem::zend = 4.0f;

#define offset1 glm::vec3(0.1f, 2.f, 0.1f)
#define offset2 glm::vec3(0.1f, 6.f, 0.1f)
/////////////////////Particle///////////////////////////////////
particle::particle(){

	}

particle::particle(glm::vec3 position){
		
		mass = 27.f;//27.f;//19.683f;=(width*radius*2)^3 * density/particle num
		force = glm::vec3(0.f);

		pos = position;
		vel = glm::vec3(0.0);

		rest_density = 1000.f;
		actual_density = rest_density;

		viscosity_coef = 100.f;//10.f
		gas_constant = 30.f;//3.f
		
		temperature = 500;

		color_interface = 1.f;
		color_surface = 1.f;

	}


//////////////////////ParticleSystem/////////////////////////////////
particleSystem::particleSystem(){

}

particleSystem::particleSystem(int number){

	initParticles(number);
	initSphere();
	
	int gx = floor(xend/SMOOTH_CORE_RADIUS) + 1;
	int gy = floor(yend/SMOOTH_CORE_RADIUS) + 1;
	int gz = floor(zend/SMOOTH_CORE_RADIUS) + 1;

	gridcells.resize(gx, gy, gz, particles);
}

void particleSystem::initParticles(int number){

	float stepsize = 2.f * radius;
	int total = number * number * number;
	
	particles.resize(total*2);
	int id;


	for(int x = 0; x < number; x++)
		for(int y = 0; y < number; y++)
			for(int z = 0; z < number; z++){
				id = x*number*number + y*number + z;
				
				particle p1(glm::vec3(x, y, z) * stepsize + offset1);
				p1.id = id;
				p1.vel.x = 5.f;
				particles[id] = p1;
				
				id += total;
				particle p2(glm::vec3(x, y, z) * stepsize + offset2);
				p2.id = id;
				particles[id] = p2;

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

void particleSystem::initCube(){
	m_positions.resize(24);
	m_colors.resize(24);
	m_normals.resize(24);
	m_indices.resize(36);

	//front
	m_positions[0] = glm::vec3( -radius,  radius,  radius );
	m_positions[1] = glm::vec3(  radius,  radius,  radius );
	m_positions[2] = glm::vec3(  radius, -radius,  radius );
	m_positions[3] = glm::vec3( -radius, -radius,  radius );
	//back
	m_positions[4] = glm::vec3( -radius,  radius, -radius );
	m_positions[5] = glm::vec3(  radius,  radius, -radius );
	m_positions[6] = glm::vec3(  radius, -radius, -radius );
	m_positions[7] = glm::vec3( -radius, -radius, -radius );
	//top
	m_positions[ 8] = glm::vec3( -radius, radius, -radius );
	m_positions[ 9] = glm::vec3(  radius, radius, -radius );
	m_positions[10] = glm::vec3(  radius, radius,  radius );
	m_positions[11] = glm::vec3( -radius, radius,  radius );
	//bottom
	m_positions[12] = glm::vec3( -radius, -radius, -radius );
	m_positions[13] = glm::vec3(  radius, -radius, -radius );
	m_positions[14] = glm::vec3(  radius, -radius,  radius );
	m_positions[15] = glm::vec3( -radius, -radius,  radius );
	//left
	m_positions[16] = glm::vec3( -radius,  radius, -radius );
	m_positions[17] = glm::vec3( -radius,  radius,  radius );
	m_positions[18] = glm::vec3( -radius, -radius,  radius );
	m_positions[19] = glm::vec3( -radius, -radius, -radius );
	//right
	m_positions[20] = glm::vec3( radius,  radius,  radius );
	m_positions[21] = glm::vec3( radius,  radius, -radius );
	m_positions[22] = glm::vec3( radius, -radius, -radius );
	m_positions[23] = glm::vec3( radius, -radius,  radius );

	//front
	m_normals[0] = glm::vec3(  0.f,  0.f,  1.f );
	m_normals[1] = glm::vec3(  0.f,  0.f,  1.f );
	m_normals[2] = glm::vec3(  0.f,  0.f,  1.f );
	m_normals[3] = glm::vec3(  0.f,  0.f,  1.f );
	//back
	m_normals[4] = glm::vec3(  0.f,  0.f, -1.f );
	m_normals[5] = glm::vec3(  0.f,  0.f, -1.f );
	m_normals[6] = glm::vec3(  0.f,  0.f, -1.f );
	m_normals[7] = glm::vec3(  0.f,  0.f, -1.f );
	//top
	m_normals[ 8] = glm::vec3(  0.f, 1.f,  0.f );
	m_normals[ 9] = glm::vec3(  0.f, 1.f,  0.f );
	m_normals[10] = glm::vec3(  0.f, 1.f,  0.f );
	m_normals[11] = glm::vec3(  0.f, 1.f,  0.f );
	//bottom
	m_normals[12] = glm::vec3(  0.f, -1.f,  0.f );
	m_normals[13] = glm::vec3(  0.f, -1.f,  0.f );
	m_normals[14] = glm::vec3(  0.f, -1.f,  0.f );
	m_normals[15] = glm::vec3(  0.f, -1.f,  0.f );
	//left
	m_normals[16] = glm::vec3( -1.f,  0.f, 0.f );
	m_normals[17] = glm::vec3( -1.f,  0.f, 0.f );
	m_normals[18] = glm::vec3( -1.f,  0.f, 0.f );
	m_normals[19] = glm::vec3( -1.f,  0.f, 0.f );
	//right
	m_normals[20] = glm::vec3( 1.f,  0.f,  0.f );
	m_normals[21] = glm::vec3( 1.f,  0.f, 0.f );
	m_normals[22] = glm::vec3( 1.f,  0.f,  0.f );
	m_normals[23] = glm::vec3( 1.f,  0.f,  0.f );


		//front
	m_indices[0 ]= 0;
	m_indices[1 ]= 1;
	m_indices[2 ]= 2;
	m_indices[3 ]= 3;
	m_indices[4 ]= 0;
	m_indices[5 ]= 2;
	//back
	m_indices[6 ]= 4;
	m_indices[7 ]= 5;
	m_indices[8 ]= 6;
	m_indices[9 ]= 7;
	m_indices[10]= 4;
	m_indices[11]= 6;
	//left
	m_indices[12]= 8;
	m_indices[13]= 9;
	m_indices[14]= 10;
	m_indices[15]= 11;
	m_indices[16]= 8;
	m_indices[17]= 10;
	//right
	m_indices[18]= 12;
	m_indices[19]= 13;
	m_indices[20]= 14;
	m_indices[21]= 15;
	m_indices[22]= 12;
	m_indices[23]= 14;
	//top
	m_indices[24]= 16;
	m_indices[25]= 17;
	m_indices[26]= 18;
	m_indices[27]= 19;
	m_indices[28]= 16;
	m_indices[29]= 18;
	//bottom
	m_indices[30]= 20;
	m_indices[31]= 21;
	m_indices[32]= 22;
	m_indices[33]= 23;
	m_indices[34]= 20;
	m_indices[35]= 22;
	
}

//how to handle boundary case?
//how to speed up neighboring search?
//what should i do if actual_density is 0?
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
		target[i].pos = source[i].pos + source[i].vel * dt 
						+ halfdt * dt * source[i].force / source[i].actual_density;
	}

	//calculate actual density 
	gridcells.refillGrid(target);

	particleGrid nghrs;
	for (int i=0; i < target.size(); i++){
		if( !checkIfOutOfBoundry(target[i]) ){
			gridcells.getNeighbors(target, target[i], nghrs);
			
			computeDensity(nghrs, target[i]);
			assert(target[i].actual_density>0);

			target[i].pressure = target[i].gas_constant * (target[i].actual_density - target[i].rest_density);
		}
	}

	

	for (int i=0; i < target.size(); i++){
		if( !checkIfOutOfBoundry(target[i]) ){
			gridcells.getNeighbors(target, target[i], nghrs);
			computeForce(nghrs, target[i]);
		}
		
	}

	//calculate actual density 
	/*for (int i=0; i < target.size(); i++){
		target[i].actual_density = computeDensity(target, target[i]);
		target[i].pressure = target[i].gas_constant * (target[i].actual_density - target[i].rest_density);
	}

	for (int i=0; i < target.size(); i++){
		target[i].force = computeForce(target, target[i]);
	}*/

	glm::vec3 collision_normal = glm::vec3(0.f);
	for (int i=0; i < target.size(); i++){
		if( CollisionDectection(target[i], collision_normal) ){
			glm::vec3 vn = (source[i].vel * collision_normal) * collision_normal;//decompose v along normal
			glm::vec3 vt = source[i].vel - vn;
			source[i].vel = 0.9f * vt - 0.5f * vn;//flip normal direction speed
		}else{
			source[i].vel += halfdt * (target[i].force/target[i].actual_density  + source[i].force /source[i].actual_density);
			source[i].pos = target[i].pos;
			source[i].force = target[i].force;
			source[i].actual_density = target[i].actual_density;
		}
	}

	

}

bool particleSystem::checkIfOutOfBoundry(const particle& p){
	glm::vec3 pos = p.pos;
	if(pos.x < xstart || pos.x > xend || pos.y < ystart || pos.y > yend || pos.z < zstart || pos.z > zend)
		return true;
	return false;
}

bool particleSystem::CollisionDectection(const particle& p, glm::vec3& n){
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
		float t = -45.f * ( ( h*h + rLen*rLen ) / rLen - 2.f * h ) / (M_PI * pow(h, 6));
		return t * r;
	}
	return glm::vec3(0.f);
}

inline float viscosityKernel(glm::vec3 r, float h){
	float rLen = glm::length(r);

	if(rLen>0 && rLen <= h)
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
void particleSystem::computeForce(const particleGrid& ps, particle& pi){
	glm::vec3 f_pressure(0.f);
	glm::vec3 f_viscosity(0.f);
	glm::vec3 f_surfaceTension(0.f);
	glm::vec3 f_gravity = glm::vec3(0,-9.8f,0) * pi.actual_density;
	
	float massOverDensity;
	glm::vec3 r;

	
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

	
	float Cs_normal_len = glm::length(Cs_normal);
	if(Cs_normal_len > surfaceThreshold){
		float curvature =  - Cs_Laplacian / Cs_normal_len;
		f_surfaceTension = tension_coeff * curvature * Cs_normal;
	}
	
	pi.force = f_pressure + f_viscosity + f_surfaceTension + f_gravity;
	//pi.force = f_pressure +  f_gravity;
}

void particleSystem::computeDensity(const particleGrid& ps, particle& pi){

	//resrt
	float rho = 0.f;

	//accumulate
	for (int j=0; j< ps.size(); j++)
	{
		rho  += ps[j].mass * poly6Kernel(pi.pos - ps[j].pos, SMOOTH_CORE_RADIUS);
	}

	pi.actual_density = rho;

		
}

///////////////////////////////draw related//////////////////////////////////////
void particleSystem::Draw(const VBO& vbos){
	
	LeapfrogIntegrate(0.01f);
	
	

	for (std::vector<particle>::iterator it = particles.begin() ; it != particles.end(); ++it){

		for(int i=0; i<= m_positions.size(); i++){
				m_positions[i] += (it->pos);
				if(it->id < particles.size() * 0.5)
					m_colors[i] = glm::vec3(0.2f,0.5f, 1.f);
				else
					m_colors[i] = glm::vec3(0.67f, 0.4f, 0.92f);
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

		for(int i=0; i<= m_positions.size(); i++){
					m_positions[i] -= (it->pos);
					m_colors[i] = glm::vec3(0.2f,0.5f, 1.f);
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
	for (std::vector<particle>::iterator it = particles.begin() ; it != particles.end(); ++it)
	{
		center = it->pos;
		io_out<<center.x<<" "<<center.y<<" "<<center.z<<" ";
	}
	io_out<<std::endl;
	i_frame++;
}
///////////////////////////////Grid//////////////////////////////////////
void particleSystem::Grid::resize(int x, int y, int z, const particleGrid& ps){
	dim = glm::vec3(x,y,z);
	GridData.resize(dim.x * dim.y * dim.z);
	refillGrid(ps);
}

void particleSystem::Grid::refillGrid(const particleGrid& ps){
	for(int i=0; i < GridData.size(); i++){
		GridData[i].clear();
		GridData[i].reserve(128);
	}

	for(int i=0; i < ps.size(); i++){

		pushParticle(ps[i]);
	}
}

void particleSystem::Grid::pushParticle(const particle& p){
	particle np = p;
	int index = positionToVecIndex(np.pos);
	/*if(p.id == 122){
		std::cout<<p.pos.x<<","<<p.pos.y<<","<<p.pos.z<<"--->";
		std::cout<<index<<std::endl;
	}*/
	if(index == -1)
		return;
	GridData[index].push_back(np.id);
}
void particleSystem::Grid::getNeighbors(const particleGrid& ps, const particle& p, particleGrid& des){
	

	des.clear();
	des.reserve(1028);
	std::vector<int> pgTemp;
	glm::vec3 gridIndex = positionToGridIndex(p.pos);
	glm::vec3 currGridIndex;
	int vecIndex;
	int pgSize;

	for(int x=-1; x<=1; x++)
		for(int y=-1; y<=1; y++)
			for(int z=-1; z<=1; z++){
				currGridIndex = gridIndex + glm::vec3(x,y,z);
				if(IfWithinBoundry(currGridIndex)){
					vecIndex = gridIndexToVecIndex(currGridIndex);
					pgTemp = GridData[vecIndex];
					for(int i=0; i < pgTemp.size(); i++)
						des.push_back(ps[pgTemp[i]]);
				}
			}
	assert(des.size()>0);

}

glm::vec3 particleSystem::Grid::positionToGridIndex(glm::vec3 p){
	glm::vec3 gridIndex(0);

	gridIndex.x = floor(p.x/SMOOTH_CORE_RADIUS);
	gridIndex.y = floor(p.y/SMOOTH_CORE_RADIUS);
	gridIndex.z = floor(p.z/SMOOTH_CORE_RADIUS);

	return gridIndex;
}


int particleSystem::Grid::gridIndexToVecIndex(glm::vec3 index){

	return dim.x * dim.y * index.z + dim.x * index.y + index.x;

}

int particleSystem::Grid::positionToVecIndex(glm::vec3 p){
	glm::vec3 gridIndex = positionToGridIndex(p);
	if( IfWithinBoundry( gridIndex ) )
		return gridIndexToVecIndex(gridIndex);
	else
		return -1;
}

bool particleSystem::Grid::IfWithinBoundry(glm::vec3 gridIndex){
	if(gridIndex.x < 0 || gridIndex.x >= dim.x)
		return false;
	if(gridIndex.y < 0 || gridIndex.y >= dim.y)
		return false;
	if(gridIndex.z < 0 || gridIndex.z >= dim.z)
		return false;
	
	return true;
}
