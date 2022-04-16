/*
  FLUIDS v.1 - SPH Fluid Simulator for CPU and GPU
  Copyright (C) 2008. Rama Hoetzlein, http://www.rchoetzlein.com

  ZLib license
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
	 claim that you wrote the original software. If you use this software
	 in a product, an acknowledgment in the product documentation would be
	 appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
	 misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
*/
#include "gl_helper.h"
#include <gl/glut.h>
#include <glm/gtx/norm.hpp>
#include <glm/gtx/string_cast.hpp>

#include "common_defs.h"
#include "fluid_system.h"

FluidSystem::FluidSystem() {
}

void FluidSystem::Draw(float* view_mat) {
	glEnable(GL_NORMALIZE);
	glLoadMatrixf(view_mat);
	for (auto& f : fluidPs) {
		glm::dvec3 scaledPos = f->pos;
		scaledPos /= SPH_RADIUS;
		glPushMatrix();
		glTranslatef(scaledPos.x, scaledPos.y, scaledPos.z);
		glScalef(SPH_RADIUS, SPH_RADIUS, SPH_RADIUS);
		glColor4f(RED(f->clr), GRN(f->clr), BLUE(f->clr), ALPH(f->clr));
		drawSphere();
		glPopMatrix();
	}
}

void FluidSystem::SPH_DrawDomain()
{
	Vector3DF min, max;
	min = m_Vec[SPH_VOLMIN];
	max = m_Vec[SPH_VOLMAX];
	min.z += 0.5;

	glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINES);
	glVertex3f(min.x, min.y, min.z);	glVertex3f(max.x, min.y, min.z);
	glVertex3f(min.x, max.y, min.z);	glVertex3f(max.x, max.y, min.z);
	glVertex3f(min.x, min.y, min.z);	glVertex3f(min.x, max.y, min.z);
	glVertex3f(max.x, min.y, min.z);	glVertex3f(max.x, max.y, min.z);
	glEnd();
}

	glVertex3f(SPH_VOLMIN.x, SPH_VOLMAX.y, SPH_VOLMIN.z);
	glVertex3f(SPH_VOLMAX.x, SPH_VOLMAX.y, SPH_VOLMIN.z);

void FluidSystem::Initialize ( int total )
{
	FreeBuffers ();
	AddBuffer(BFLUID, sizeof(Fluid), total);
	//for (int i = 0; i < total; ++i) {
	//	fluidPs.push_back(std::make_unique<Fluid>());
	//}
		
	SPH_Setup ();
	Reset ( total );	
}

void FluidSystem::Reset ( int nmax )
{
	ResetBuffer ( 0, nmax );

	m_DT = 0.003; //  0.001;			// .001 = for point grav

	// Reset parameters
	m_Param [ MAX_FRAC ] = 1.0;
	m_Param [ POINT_GRAV ] = 0.0;
	m_Param [ PLANE_GRAV ] = 1.0;

	m_Param [ BOUND_ZMIN_SLOPE ] = 0.0;
	m_Param [ FORCE_XMAX_SIN ] = 0.0;
	m_Param [ FORCE_XMIN_SIN ] = 0.0;	
	m_Param [ SPH_INTSTIFF ] = 1.00;
	m_Param [ SPH_VISC ] = 0.2;
	m_Param [ SPH_INTSTIFF ] = 0.50;
	m_Param [ SPH_EXTSTIFF ] = 20000;
	m_Param [ SPH_SMOOTHRADIUS ] = 0.01;
	
	m_Vec [ POINT_GRAV_POS ].Set ( 0, 0, 50 );
	m_Vec [ PLANE_GRAV_DIR ].Set ( 0, 0, -9.8 );
}

int FluidSystem::AddPoint ()
{
	xref ndx;	
	Fluid* f = (Fluid*) AddElem ( 0, ndx );	
	f->sph_force.Set(0,0,0);
	f->vel.Set(0,0,0);
	f->vel_eval.Set(0,0,0);
	f->next = 0x0;
	f->pressure = 0;
	f->density = 0;
	return ndx;
}

int FluidSystem::AddPointReuse ()
{
	xref ndx;	
	Fluid* f;
	if ( NumPoints() <= mBuf[0].max-2 )
		f = (Fluid*) AddElem ( 0, ndx );
	else
		f = (Fluid*) RandomElem ( 0, ndx );

	f->sph_force.Set(0,0,0);
	f->vel.Set(0,0,0);
	f->vel_eval.Set(0,0,0);
	f->next = 0x0;
	f->pressure = 0;
	f->density = 0;
	return ndx;
}

void FluidSystem::Run ()
{	
	Grid_InsertParticles ();
	SPH_ComputePressureGrid ();
	SPH_ComputeForceGridNC ();
	Advance();
}

void FluidSystem::Advance ()
{
	char *dat1, *dat1_end;
	Fluid* p;
	Vector3DF norm, z;
	Vector3DF dir, accel;
	Vector3DF vnext;
	Vector3DF min, max;
	double adj;
	float SL, SL2, ss, radius;
	float stiff, damp, speed, diff; 
	SL = m_Param[SPH_LIMIT];
	SL2 = SL*SL;
	
	stiff = m_Param[SPH_EXTSTIFF];
	damp = m_Param[SPH_EXTDAMP];
	radius = m_Param[SPH_PRADIUS];
	min = m_Vec[SPH_VOLMIN];
	max = m_Vec[SPH_VOLMAX];
	ss = m_Param[SPH_SIMSCALE];

	dat1_end = mBuf[0].data + NumPoints()*mBuf[0].stride;
	for ( dat1 = mBuf[0].data; dat1 < dat1_end; dat1 += mBuf[0].stride ) {
		p = (Fluid*) dat1;		

		// Compute Acceleration		
		accel = p->sph_force;
		accel *= m_Param[SPH_PMASS];

		// Velocity limiting 
		speed = accel.x*accel.x + accel.y*accel.y + accel.z*accel.z;
		if ( speed > SL2 ) {
			accel *= SL / sqrt(speed);
		}		
	
		// Boundary Conditions
		// Z-axis walls (floor)
		diff = 2 * radius - ( p->pos.z - min.z - (p->pos.x - m_Vec[SPH_VOLMIN].x) * m_Param[BOUND_ZMIN_SLOPE] )*ss;
		if (diff > EPSILON) {
			norm.Set(-m_Param[BOUND_ZMIN_SLOPE], 0, 1.0 - m_Param[BOUND_ZMIN_SLOPE]);
			adj = stiff * diff - damp * norm.Dot(p->vel_eval);
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}

		diff = 2 * radius - ( max.z - p->pos.z )*ss;
		if (diff > EPSILON) {
			norm.Set ( 0, 0, -1 );
			adj = stiff * diff - damp * norm.Dot ( p->vel_eval );
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}
		
		// X-axis walls
		diff = 2 * radius - ( p->pos.x - min.x + (sin(m_Time*10.0)-1+(p->pos.y*0.025)*0.25) * m_Param[FORCE_XMIN_SIN] )*ss;	
		//diff = 2 * radius - ( p->pos.x - min.x + (sin(m_Time*10.0)-1) * m_Param[FORCE_XMIN_SIN] )*ss;	
		if (diff > EPSILON ) {
			norm.Set ( 1.0, 0, 0 );
			adj = (m_Param[ FORCE_XMIN_SIN ] + 1) * stiff * diff - damp * norm.Dot ( p->vel_eval ) ;
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;					
		}

int FluidSystem::GetGridIndex(const glm::ivec3 &gridPos) {
	return gridPos.z * gridSpaceDiag.y * gridSpaceDiag.x + gridPos.y * gridSpaceDiag.x + gridPos.x;
}

		// Y-axis walls
		diff = 2 * radius - ( p->pos.y - min.y )*ss;			
		if (diff > EPSILON) {
			norm.Set ( 0, 1, 0 );
			adj = stiff * diff - damp * norm.Dot ( p->vel_eval );
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}
		diff = 2 * radius - ( max.y - p->pos.y )*ss;
		if (diff > EPSILON) {
			norm.Set ( 0, -1, 0 );
			adj = stiff * diff - damp * norm.Dot ( p->vel_eval );
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}

		// Plane gravity
		if ( m_Param[PLANE_GRAV] > 0) 
			accel += m_Vec[PLANE_GRAV_DIR];

		// Point gravity
		if ( m_Param[POINT_GRAV] > 0 ) {
			norm.x = ( p->pos.x - m_Vec[POINT_GRAV_POS].x );
			norm.y = ( p->pos.y - m_Vec[POINT_GRAV_POS].y );
			norm.z = ( p->pos.z - m_Vec[POINT_GRAV_POS].z );
			norm.Normalize ();
			norm *= m_Param[POINT_GRAV];
			accel -= norm;
		}

		p->predictPos = p->pos + (p->vel * m_DT);


		// Perform collision detection and response
		if (p->predictPos.y < scaledMin.y) { p->vel.y = 0.0; p->predictPos.y = scaledMin.y + 0.001; }
		if (p->predictPos.y > scaledMax.y) { p->vel.y = 0.0; p->predictPos.y = scaledMax.y - 0.001; }

		if (p->predictPos.x < scaledMin.x) { p->vel.x = 0.0; p->predictPos.x = scaledMin.x + 0.001; }
		if (p->predictPos.x > scaledMax.x) { p->vel.x = 0.0; p->predictPos.x = scaledMax.x - 0.001; }

		if (p->predictPos.z < scaledMin.z) { p->vel.z = 0.0; p->predictPos.z = scaledMin.z + 0.001; }
		if (p->predictPos.z > scaledMax.z) { p->vel.z = 0.0; p->predictPos.z = scaledMax.z - 0.001; }
	}
}

void FluidSystem::FindNeighbors() {
	for (int i = 0; i < totalGridCells; ++i) {
		grid.at(i).clear();
	}

	//equivalinet of insertgrid / update grid finding the postns within the grid
	for (int i = 0; i < fluidPs.size(); ++i) {
		glm::ivec3 gridPos = GetGridPos(fluidPs.at(i)->predictPos);
		int gIndex = GetGridIndex(gridPos);

		// this if shouldn't be necessary?
		if (0 <= gIndex && gIndex < grid.size()) {
			grid.at(gIndex).push_back(i); // maybe set a limit? (see MAX_NEIGHBOR)
		}
	}

	// equiv. Finding the Neighbors
	for (int i = 0; i < fluidPs.size(); ++i) {
		neighbor_loop:
		neighbors.at(i).clear(); // clear neighbors from prev.

		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		glm::ivec3 gridPos = GetGridPos(p->predictPos);

		int SEARCH_SIZE = 1;
		// 2x2 neighborhood.
		for (int x = -SEARCH_SIZE; x <= SEARCH_SIZE; x++) {
		for (int y = -SEARCH_SIZE; y <= SEARCH_SIZE; y++) {
		for (int z = -SEARCH_SIZE; z <= SEARCH_SIZE; z++) {
			glm::ivec3 n = gridPos + glm::ivec3(x, y, z);
			if (0 <= n.x && n.x < gridSpaceDiag.x &&
				0 <= n.y && n.y < gridSpaceDiag.y &&
				0 <= n.z && n.z < gridSpaceDiag.z) {
				int gIndex = GetGridIndex(n);
				for (int pIndex : grid.at(gIndex)) { // each 
					std::unique_ptr<Fluid>& pcurr = fluidPs.at(pIndex);
					double lenR = glm::length(p->predictPos - pcurr->predictPos);
					if (lenR <= SPH_RADIUS) {
						neighbors.at(i).push_back(pIndex);
						if (neighbors.at(i).size() >= MAX_NEIGHBOR) { goto neighbor_loop; }
					}
				}
			}
		}}}
	}
}

void FluidSystem::ComputeDensity() {
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		p->density = 0.0;
		for (int j : neighbors.at(i)) { // for each neighbor
			p->density += PolyKernel(glm::length(p->predictPos - fluidPs.at(j)->predictPos));
		}
	}
}

void FluidSystem::ComputeLambda() {
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		double sumGradients = 0.0;
		glm::dvec3 pGrad = glm::dvec3(0.0);
		for (int j : neighbors.at(i)) { // for each neighbor
			std::unique_ptr<Fluid>& pcurr = fluidPs.at(j);

			// Spiky Kernel - modifies r by ref
			glm::dvec3 r = (p->predictPos - pcurr->predictPos);
			SpikyKernel(r);
			r /= REST_DENSITY;
			// End Spiky Kernel
			sumGradients += glm::length2(r);
			pGrad += r; // -= r; ?? - i think += b/c -45
		}
		sumGradients += glm::length2(pGrad);
		double constraint = p->density / REST_DENSITY - 1.0; // real scale constraint
		p->lambda = -constraint / (sumGradients + RELAXATION); // maybe + 500 or so
	}
}

void FluidSystem::ComputeCorrections() {
	double polyDen = PolyKernel(0.2 * SPH_RADIUS);
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		p->deltaPos = glm::dvec3(0.0);
		for (int j : neighbors.at(i)) { // for each neighbor
			std::unique_ptr<Fluid>& pcurr = fluidPs.at(j);
			//---------Calculate SCORR-----
			double frac = PolyKernel(glm::length(p->predictPos - pcurr->predictPos)) / polyDen;
			double sCorr = -K_CORR * frac * frac * frac * frac;
			//------------End SCORR calculation-------

			// Spiky Kernel - modifies r by ref
			glm::dvec3 r = (p->predictPos - pcurr->predictPos);
			SpikyKernel(r);
			r /= REST_DENSITY;
			// End Spiky Kernel

			//p->deltaPos += r * (p->lambda + pcurr->lambda
			p->deltaPos += r * (p->lambda + pcurr->lambda + sCorr);
		}
	}
}

void FluidSystem::ApplyCorrections() {
	for (std::unique_ptr<Fluid>& p : fluidPs) {
		p->predictPos += p->deltaPos;
	}
}

void FluidSystem::Advance() {
	//update all velocities
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		//std::cout << "for particle: " << i << std::endl;
		//std::cout << "p pos: " << glm::to_string(p->pos) << std::endl;
		//std::cout << "p pred: " << glm::to_string(p->predictPos) << std::endl;
		//std::cout << "p delP: " << glm::to_string(p->deltaPos) << std::endl;
		//std::cout << "p vel: " << glm::to_string(p->vel) << std::endl;
		//std::cout << "p dens: " << p->density << std::endl;
		//std::cout << "p lamb: " << p->lambda << std::endl << std::endl;
		p->vel = (p->predictPos - p->pos) / m_DT;

		// vorticity confinement here?
		p->pos = p->predictPos;
	}

	// VISCOSITY
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		glm::dvec3 acc(0.0, 0.0, 0.0);
		for (int j : neighbors.at(i)) {
			std::unique_ptr<Fluid>& pcurr = fluidPs.at(j);
			acc += (pcurr->vel - p->vel) * PolyKernel(glm::length(p->predictPos - pcurr->predictPos));
		}
		p->vel_tmp = acc;
	}
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		p->vel += VISC_CONST * p->vel_tmp * m_DT;
	}
	// END VISCOSITY
}