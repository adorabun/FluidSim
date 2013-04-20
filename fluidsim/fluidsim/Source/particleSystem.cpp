#include "particleSystem.h"
#include <cmath>
#include <fstream>
#include <iostream>

//#define EPSILON 0.000001f
#define SMOOTH_CORE_RADIUS 0.5f //three times the average distance


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

#define offset1 glm::vec3(2.2f, 2.f, 3.5f)
#define offset2 glm::vec3(2.2f, 5.5f, 3.5f)
/////////////////////Particle///////////////////////////////////
particle::particle(){

	}

particle::particle(glm::vec3 position, float rho){
		
		rest_density = rho;
		actual_density = rest_density;

		mass = rho * 1.333333f * M_PI * pow(particleSystem::radius,3);//27.f;//19.683f;=(width*radius*2)^3 * density/particle num
		force = glm::vec3(0.f);

		pos = position;
		vel = glm::vec3(0.0);

		

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

	container = new Mesh();
	container->readOBJ("obj/capsule_half.obj");
	container->scale(2.5f);

	frameCount = 0;

	initParticles(number);
	initSphere();
	
	int gx = floor(10.0/SMOOTH_CORE_RADIUS) + 1;
	int gy = floor(30.0/SMOOTH_CORE_RADIUS) + 1;
	int gz = floor(10.0/SMOOTH_CORE_RADIUS) + 1;

	//gridcells.resize(gx, gy, gz, particles);
	mygrid.dim = glm::vec3(gx,gy,gz);
}
//void particleSystem::initParticles(int number){
//
//	float stepsize = 2.f * radius;
//	int total = number * number * number;
//	
//	particles.resize(total);
//	
//	int id;
//
//
//	for(int x = 0; x < number; x++)
//		for(int y = 0; y < number; y++)
//			for(int z = 0; z < number; z++){
//				id = x*number*number + y*number + z;
//				
//				particle* p1 = new particle(glm::vec3(x, y, z) * stepsize + offset1);
//				p1->id = id;
//				p1->vel.y = -5.f;
//				particles[id] = p1;
//
//			}
//
//}

//void particleSystem::initParticles(int number){
//
//	float stepsize = 2.f * radius;
//	int total = number * number * number;
//	
//	particles.resize(total);
//	
//	int id;
//
//
//	for(int x = 0; x < number; x++)
//		for(int y = 0; y < number; y++)
//			for(int z = 0; z < number; z++){
//				id = x*number*number + y*number + z;
//				
//				particle* p1 = new particle(glm::vec3(x, y, z) * stepsize + offset1);
//				p1->id = id;
//				p1->vel.y = -5.f;
//				particles[id] = p1;
//
//			}
//
//}

void particleSystem::initParticles(int number){

	float stepsize = 2.f * radius;
	int total = number * number * number;
	
	particles.resize(total*2);
	
	int id;


	for(int x = 0; x < number; x++)
		for(int y = 0; y < number; y++)
			for(int z = 0; z < number; z++){
				id = x*number*number + y*number + z;
				
				particle p1(glm::vec3(x, y, z) * stepsize + offset1, 500.f);
				p1.id = id;
				particles[id] = p1;
				
				id += total;
				particle p2(glm::vec3(x, y, z) * stepsize + offset2, 1000.f);
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
			m_positions[index] = glm::vec3(x, y, z) * radius * 0.8f;;
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


void particleSystem::LeapfrogIntegrate(float dt){
	float halfdt = 0.5f * dt;
	std::vector<particle> target = particles;// target is a copy!
	std::vector<particle>& source = particles;//source is a ptr!
	glm::vec3 collision_normal = glm::vec3(0.f);

	
	//double time0 = glfwGetTime();
	for (int i=0; i < target.size(); i++){
		target[i].pos = source[i].pos + source[i].vel * dt 
						+ halfdt * dt * source[i].force / source[i].actual_density;
		if(container->lineIntersect(source[i].pos, target[i].pos, collision_normal) ){
		//if( CollisionDectection(target[i], collision_normal) ){
			glm::vec3 vn = (source[i].vel * collision_normal) * collision_normal;//decompose v along normal
			glm::vec3 vt = source[i].vel - vn;
			target[i].vel = 0.9f * vt - 0.8f * vn;//flip normal direction speed
			target[i].pos = source[i].pos;
		}

		mygrid.pushParticle(target[i], frameCount);

	}

	//double time1 = glfwGetTime();

	
	for (int i=0; i < target.size(); i++){
		mygrid.getNeighbors(target[i]);
	}
	
	//double time2= glfwGetTime();

	//calculate actual density 
	for (int i=0; i < target.size(); i++){
		//mygrid.getNeighbors(target[i]);
			
		computeDensity(target[i]);
		

		target[i].pressure = target[i].gas_constant * (target[i].actual_density - target[i].rest_density);
	}

	
	//double time3= glfwGetTime();
	for (int i=0; i < target.size(); i++){
		
		computeForce(target[i]);	
	}
	//double time4= glfwGetTime();

	//calculate actual density 
	/*for (int i=0; i < target.size(); i++){
		computeDensity(target, target[i]);
		target[i]->pressure = target[i]->gas_constant * (target[i]->actual_density - target[i]->rest_density);
	}

	for (int i=0; i < target.size(); i++){
		computeForce(target, target[i]);
	}*/

	
	for (int i=0; i < target.size(); i++){
		source[i].pos = target[i].pos;
		source[i].vel = target[i].vel + halfdt * (target[i].force/target[i].actual_density  + source[i].force /source[i].actual_density);
		source[i].force = target[i].force;
		source[i].actual_density = target[i].actual_density;	
	}

	//double time5 = glfwGetTime();

	/*std::cout<<"============"<<frameCount<<std::endl
			 << "push   Particle = 2-1=" << time1-time0<<std::endl
			 << "get    Neighbor = 2-1=" << time2-time1<<std::endl
			 << "compute Density = 3-2=" << time3-time2<<std::endl
			 << "compute   Force = 4-3=" << time4-time3<<std::endl
			 << "copy       back = 5-4=" << time5-time4<<std::endl;*/
	frameCount++;
}

bool particleSystem::checkIfOutOfBoundry(const particle& p){
	glm::vec3 pos = p.pos;
	if(pos.x < xstart || pos.x > xend || pos.y < ystart || pos.y > yend || pos.z < zstart || pos.z > zend)
		return true;
	return false;
}

bool particleSystem::CollisionDectection(particle& p, glm::vec3& n){
	n = glm::vec3(0.f);
	glm::vec3 pos = p.pos;

	if(pos.x < xstart){
		n.x = 1.f;
		//p.pos.x = xstart + EPSILON;
	}
	else if(pos.x > xend){
		n.x = -1.f;
		//p.pos.x = xend - EPSILON;
	}
	
	if(pos.y < ystart){
		n.y = 1.f;
		//p.pos.y = ystart + EPSILON;
	}
	else if(pos.y > yend){
		n.y = -1.f;
		//p.pos.y = yend - EPSILON;
	}

	if(pos.z < zstart){
		n.z = 1.f;
		//p.pos.z = zstart + EPSILON;
	}
	else if(pos.z > zend){
		n.z = -1.f;
		//p.pos.z = zend - EPSILON;
	}

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
		float t =  - 945.f  * pow(h*h - rLen*rLen, 2) / (32.f * M_PI * pow(h,9) );
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

	if(rLen>0 && rLen <= h)
		return 15.f * ( - rLen*rLen*rLen/(2.f*h*h*h) + rLen*rLen/(h*h) + h/(2.f*rLen) - 1 ) / ( 2 * M_PI * pow( h, 3 ) ) ;
	return 0.f;
}

inline float viscosityKernelLaplacian(glm::vec3 r, float h){
	float rLen = glm::length(r);
	if(rLen <= h)
		return 45.f * ( h - rLen ) / ( M_PI * pow(h, 6) ) ;
	return 0.f;
}

///////////////////////////////computation//////////////////////////////////////
void particleSystem::computeForce(particle& pi){

	glm::vec3 f_pressure(0.f);
	glm::vec3 f_viscosity(0.f);
	glm::vec3 f_surfaceTension(0.f);
	glm::vec3 f_gravity = glm::vec3(0,-9.8f * 3,0) * pi.actual_density;
	
	float massOverDensity;
	glm::vec3 r;

	
	glm::vec3 Cs_normal = glm::vec3(0.f);
	float Cs_Laplacian = 0.f;

	for (int j=0; j< pi.ngbrs.size(); j++){
		massOverDensity = pi.ngbrs[j]->mass / pi.ngbrs[j]->actual_density;

		r = pi.pos - pi.ngbrs[j]->pos;

		
		f_pressure -=  massOverDensity * ( pi.pressure + pi.ngbrs[j]->pressure ) * 0.5f * spikyKernelGradient(r, SMOOTH_CORE_RADIUS);

		f_viscosity  += ( pi.viscosity_coef + pi.ngbrs[j]->viscosity_coef) * 0.5f *
			massOverDensity * ( pi.ngbrs[j]->vel - pi.vel ) * viscosityKernelLaplacian(r, SMOOTH_CORE_RADIUS);

		Cs_normal += massOverDensity * poly6KernelGradient(r, SMOOTH_CORE_RADIUS);
		Cs_Laplacian += massOverDensity * poly6KernelLaplacian(r, SMOOTH_CORE_RADIUS);
		
	}

	float Cs_normal_len = glm::length(Cs_normal);
	if(Cs_normal_len > surfaceThreshold){
		float curvature =  - Cs_Laplacian / Cs_normal_len;
		f_surfaceTension = tension_coeff * curvature * Cs_normal;
	}
	
	pi.force = f_pressure + f_viscosity + f_surfaceTension + f_gravity;
	//pi.force = f_pressure +  f_gravity;
}

void particleSystem::computeDensity(particle& pi){
	
	//resrt
	float rho = 0.f;

	bool containself = false;
	//accumulate
	for (int j=0; j< pi.ngbrs.size(); j++)
	{
		if(pi.ngbrs[j]->id == pi.id)
			containself = true;
		rho  += pi.ngbrs[j]->mass * poly6Kernel(pi.pos - pi.ngbrs[j]->pos, SMOOTH_CORE_RADIUS);
	}

	/*if(!containself)
		std::cout << "framecount =" <<frameCount <<" id= "<<pi->id<<std::endl;*/

	//if( rho <= 0.f){
	//	std::cout <<"current id = "<< pi->id << std::endl;
	//	for (int j=0; j< pi->ngbrs.size(); j++)
	//		std::cout <<"count=" << j << " id=" << (pi->ngbrs[j]->id) << " mass =" <<pi->ngbrs[j]->mass 
	//		<<" distance= " << glm::length(pi->pos - pi->ngbrs[j]->pos)
	//				<<" kernel= " << poly6Kernel(pi->pos - pi->ngbrs[j]->pos, SMOOTH_CORE_RADIUS)<<std::endl;

	//}
	
	pi.actual_density = rho;

	assert(rho > 0);	
}

///////////////////////////////draw related//////////////////////////////////////
void particleSystem::Draw(const VBO& vbos){

	LeapfrogIntegrate(0.01f);
	
	for (int id=0; id<particles.size(); id++ ){

		for(int i=0; i< m_positions.size(); i++){
				m_positions[i] += (particles[id].pos);
				/*if(id==37)
					m_colors[i] = glm::vec3(1.f,0,0);
				else */if(id> particles.size() * 0.5)
					m_colors[i] = glm::vec3(0.2f,0.5f, 1.f);//blue
				else
					m_colors[i] = glm::vec3(0.89f, 0.71f, 0.21f);//oil
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

		for(int i=0; i< m_positions.size(); i++){
					m_positions[i] -= (particles[id].pos);
					m_colors[i] = glm::vec3(0.2f,0.5f, 1.f);
		}

	}
	container->draw(vbos);
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
	for (int id=0; id<particles.size(); id++ )
	{
		center = particles[id].pos;
		io_out<<center.x<<" "<<center.y<<" "<<center.z<<" ";
	}
	io_out<<std::endl;
	i_frame++;
}


//map<int,Box> for hash
//add frame count, only clear for each frame
//draw function, bind once, only manipulate
//change to vector<particle*> checked

///////////////////////////////SpaceGrid//////////////////////////////////////


void SpaceGrid::pushParticle(particle& pt, int frameID){
	
	int index = positionToVecIndex(pt.pos);

	if(index == -1)
		return;


	if(mymap[index].frameID < frameID){
		mymap[index].frameID = frameID;
		mymap[index].ps.clear();
		mymap[index].ps.reserve(64);
	}
	
	mymap[index].ps.push_back(&pt);

}
void SpaceGrid::getNeighbors(particle& pt){
	

	std::vector<particle*> temp;
	
	glm::vec3 gridIndex = positionToGridIndex(pt.pos);
	glm::vec3 currGridIndex;
	int vecIndex;

	for(int x = -1; x <= 1; x++)
		for(int y = -1; y <= 1; y++)
			for(int z = -1; z <= 1; z++){
				currGridIndex = gridIndex + glm::vec3(x,y,z);
				
				vecIndex = gridIndexToVecIndex(currGridIndex);
				
				
				if( mymap.find(vecIndex) != mymap.end() )//sanity check
					temp.insert(temp.end(), mymap[vecIndex].ps.begin(), mymap[vecIndex].ps.end() );
					
			}

	pt.ngbrs = temp;

	//std::cout<<pt->ngbrs.size()<<std::endl;

	/*bool containitself = false;
	for(int i=0; i<temp.size(); i++){
		if( temp[i]->id == pt.id){
			containitself = true;
			break;
		}
	}

	if(!containitself){
		std::cout  <<" id= "<<pt->id<<std::endl;
		if(pt.id==37){
			gridIndex = positionToGridIndex(pt.pos);
			vecIndex = gridIndexToVecIndex(gridIndex);
			std::cout<<"pos="<<pt.pos.x <<", "<<pt.pos.y<<", "<<pt.pos.z<<std::endl 
				<<"gridIndex ="<<gridIndex.x <<", "<<gridIndex.y<<", "<<gridIndex.z<<std::endl
				<<"vecIndex = "<<vecIndex<<std::endl;

		}
	}*/

	assert(pt.ngbrs.size()>0);

}

glm::vec3 SpaceGrid::positionToGridIndex(glm::vec3 p){
	glm::vec3 gridIndex(0);

	gridIndex.x = floor(p.x/SMOOTH_CORE_RADIUS);
	gridIndex.y = floor(p.y/SMOOTH_CORE_RADIUS);
	gridIndex.z = floor(p.z/SMOOTH_CORE_RADIUS);

	return gridIndex;
}


int SpaceGrid::gridIndexToVecIndex(glm::vec3 index){
	
	if( IfWithinBoundry( index) )
		return dim.x * dim.y * index.z + dim.x * index.y + index.x;
	else
		return -1;

}

int SpaceGrid::positionToVecIndex(glm::vec3 p){
	return gridIndexToVecIndex( positionToGridIndex(p) );

}

bool SpaceGrid::IfWithinBoundry(glm::vec3 gridIndex){
	if(gridIndex.x < 0 || gridIndex.x >= dim.x)
		return false;
	if(gridIndex.y < 0 || gridIndex.y >= dim.y)
		return false;
	if(gridIndex.z < 0 || gridIndex.z >= dim.z)
		return false;
	
	return true;
}