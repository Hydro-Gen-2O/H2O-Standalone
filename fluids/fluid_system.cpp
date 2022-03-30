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
#include <conio.h>

#include <gl/glut.h>
#include <glm/gtx/norm.hpp>
#include <glm/gtx/string_cast.hpp>

#include "common_defs.h"
#include "fluid_system.h"

#define EPSILON			0.00001f			//for collision detection

FluidSystem::FluidSystem()
{
}

void FluidSystem::SPH_DrawDomain()
{
	glm::vec3 min, max;
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

void FluidSystem::SPH_Setup(int n) {
	//m_Param[BOUND_ZMIN_SLOPE] = 0.0;
	//m_Param[FORCE_XMAX_SIN] = 0.0;
	//m_Param[FORCE_XMIN_SIN] = 0.0;
	m_Param[SPH_INTSTIFF] = 0.50;
	m_Param[PLANE_GRAV] = 0.0;
	m_Param[SPH_VISC] = 0.2;			// pascal-second (Pa.s) = 1 kg m^-1 s^-1  (see wikipedia page on viscosity)

	m_Param[SPH_SMOOTHRADIUS] = 0.01; // 0.01;

	m_Vec[POINT_GRAV_POS] = glm::vec3(0, 0, 50);
	m_Vec[PLANE_GRAV_DIR] = glm::vec3(0, 0, -9.8);

	m_Param[SPH_SIMSCALE] = 0.004;// 0.004;			// unit size
	m_Param[SPH_RESTDENSITY] = 600.0; //600.0;			// kg / m^3
	m_Param[SPH_PMASS] = 0.00020543;// 0.00020543;		// kg

	m_Param[SPH_LIMIT] = 200.0;			// m / s

	m_Param[SPH_PRADIUS] = 0.004;			// m
	m_Param[SPH_EXTSTIFF] = 20000;
	m_Param[SPH_EXTDAMP] = 256.0;
	// kernel computation
	m_Poly6Kern = 315.0f / (64.0f * 3.141592 * pow(m_Param[SPH_SMOOTHRADIUS], 9));	// Wpoly6 kernel (denominator part) - 2003 Muller, p.4
	m_SpikyKern = -45.0f / (3.141592 * pow(m_Param[SPH_SMOOTHRADIUS], 6));			// grad of spiky (no minus)
	m_LapKern = 45.0f / (3.141592 * pow(m_Param[SPH_SMOOTHRADIUS], 6));

	switch (n) {
	case -1:		// fluid drop
		m_Vec[SPH_VOLMIN] = glm::vec3(-10, -10, 0);
		m_Vec[SPH_VOLMAX] = glm::vec3(10, 10, 30);
		// INIT MAX/ INIT MIN governs the two opposing corners of a square drop 
		m_Vec[SPH_INITMIN] = glm::vec3(-5, -5, 10);
		m_Vec[SPH_INITMAX] = glm::vec3(5, 5, 20);
		break;
	}

	m_DT = 0.003; //  0.001;			// .001 = for point grav

	Grid_Setup(
		m_Vec[SPH_VOLMIN],
		m_Vec[SPH_VOLMAX],
		m_Param[SPH_SIMSCALE],
		m_Param[SPH_SMOOTHRADIUS] * 2.0,
		1.0);
}

void FluidSystem::SPH_CreateExample(int n, int nmax)
{
	SPH_Setup(n);
	fluidPs.clear();
	maxPoints = nmax;
	//float volume = abs(m_Vec[SPH_INITMIN].x - m_Vec[SPH_INITMAX].x) *
	//	abs(m_Vec[SPH_INITMIN].y - m_Vec[SPH_INITMAX].y) *
	//	abs(m_Vec[SPH_INITMIN].z - m_Vec[SPH_INITMAX].z);

	//TODO: weird spawning
	// volume^1/3 * 0.87 / simscale?
	float ss = pow(m_Param[SPH_PMASS] / m_Param[SPH_RESTDENSITY], 1 / 3.0) * 0.87 / m_Param[SPH_SIMSCALE];

	glm::vec3 min = m_Vec[SPH_INITMIN];
	glm::vec3 max = m_Vec[SPH_INITMAX];
	float dx = max.x - min.x;
	float dy = max.y - min.y;
	float dz = max.z - min.z;
	for (float z = max.z; z >= min.z; z -= ss) {
		for (float y = min.y; y <= max.y; y += ss) {
			for (float x = min.x; x <= max.x; x += ss) {
				int ndx;
				//if (fluidPs.size() < maxPoints - 2) {
				ndx = fluidPs.size();
				fluidPs.push_back(std::make_unique<Fluid>());
				//}
				//else {
				//	ndx = rand() % fluidPs.size();
				//}
				fluidPs.at(ndx)->pos = glm::vec3(x, y, z);
				fluidPs.at(ndx)->predictPos = glm::vec3(x, y, z);
				fluidPs.at(ndx)->clr = COLORA((x - min.x) / dx, (y - min.y) / dy, (z - min.z) / dz, 1);
			}
		}
	}
}

void FluidSystem::Run()
{
	//PBF_PredictPositions();
	//SPH_FindNeighbors(true); //?
	//for (int i = 0; i < 4; ++i) {
	//	SPH_ComputeDensity();
	//	SPH_ComputeLambda();
	//	SPH_ComputeCorrections();
	//	SPH_ApplyCorrections();
	//}
	//Advance();

	SPH_FindNeighbors(false);
	SPH_ComputeDensity();
	SPH_ComputeForceGridNC();
	AdvanceOld();
}

void FluidSystem::PBF_PredictPositions() {
	for (std::unique_ptr<Fluid>& p : fluidPs) {
		glm::vec3 deltaVel = m_Vec[PLANE_GRAV_DIR];
		//deltaVel *= m_Param[SPH_PMASS]; // m * a
		deltaVel *= m_DT; // a * dt = change in v

		// TODO TOCHANGE NOTE apply force to velocity (gravity)
		//p->vel += deltaVel;

		// limit velocity ? 
		//if (length vec > m_Param[SPH_LIMIT]) {
		//	p->vel = m_Param[SPH_LIMIT];
		//}

		glm::vec3 velDist = p->vel * m_DT;
		// /; // ?????????
		//velDist /= m_Param[SPH_SIMSCALE];
		p->predictPos = p->pos + velDist;
	}
}

void FluidSystem::SPH_FindNeighbors(bool PBF) {
	glm::vec3 d;
	Grid_InsertParticles();
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		float sum = 0.0;
		m_NC[i] = 0;
		if (PBF) {
			Grid_FindCells(p->predictPos, m_Param[SPH_SMOOTHRADIUS] / m_Param[SPH_SIMSCALE]);
		} else {
			Grid_FindCells(p->pos, m_Param[SPH_SMOOTHRADIUS] / m_Param[SPH_SIMSCALE]);
		}
		for (int cell = 0; cell < 8; cell++) {
			if (m_GridCell[cell] != -1) {
				int pndx = m_Grid[m_GridCell[cell]];
				while (pndx != -1) {
					std::unique_ptr<Fluid>& pcurr = fluidPs.at(pndx);
					if (pndx == i) { pndx = pcurr->next; continue; }
					if (PBF) {
						d = (p->predictPos - pcurr->predictPos) * m_Param[SPH_SIMSCALE];
					}
					else {
						d = (p->pos - pcurr->pos) * m_Param[SPH_SIMSCALE];
					}
					float lenR = glm::length(d);
					if (lenR < m_Param[SPH_SMOOTHRADIUS]) {
						if (m_NC[i] < MAX_NEIGHBOR) {
							m_Neighbor[i][m_NC[i]] = pndx;
							m_NDist[i][m_NC[i]] = lenR;
							m_NC[i]++;
						}
					}
					pndx = pcurr->next;
				}
			}
			m_GridCell[cell] = -1;
		}
	}
}

void FluidSystem::SPH_ComputeDensity() {
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		float sum = 0.0;
		for (int j = 0; j < m_NC[i]; j++) {
			float c = (float)m_Param[SPH_SMOOTHRADIUS] * (float)m_Param[SPH_SMOOTHRADIUS] -
				m_NDist[i][j] * m_NDist[i][j];
			sum += c * c * c;
		}
		p->density = sum * m_Param[SPH_PMASS] * m_Poly6Kern;
		// pressure only useful for traditional SPH method
		p->pressure = (p->density - m_Param[SPH_RESTDENSITY]) * m_Param[SPH_INTSTIFF];
	}
}

void FluidSystem::SPH_ComputeLambda() {
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);

		//WHICH?
		//float constraint = p->density / m_Param[SPH_RESTDENSITY] - 1.f;
		float constraint = p->density / m_Param[SPH_RESTDENSITY] - 1.f;

		p->gradient = glm::vec3(0.0, 0.0, 0.0);
		float sumGradients = 0.0f;
		for (int j = 0; j < m_NC[i]; j++) {
			std::unique_ptr<Fluid>& pcurr = fluidPs.at(m_Neighbor[i][j]);

			// Spiky Kernel
			glm::vec3 r = (p->predictPos - pcurr->predictPos) * m_Param[SPH_SIMSCALE];
			r *= m_SpikyKern;			// grad of spiky (no minus)
			r *= (m_Param[SPH_SMOOTHRADIUS] - m_NDist[i][j]) * (m_Param[SPH_SMOOTHRADIUS] - m_NDist[i][j]) /
				m_NDist[i][j];
			r /= m_Param[SPH_RESTDENSITY];
			// End Spiky Kernel

			sumGradients += glm::length2(r);
			p->gradient += r; // -= r; ?? - i think += b/c -45
		}
		sumGradients += glm::length2(p->gradient);
		p->lambda = -constraint / (sumGradients); // maybe + 500 or so

		// BELOW - not sure ifn ecessray?
		p->lambda /= m_Param[SPH_SIMSCALE];
	}
}

void FluidSystem::SPH_ComputeCorrections() {
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		p->deltaPos = glm::vec3(0.f, 0.f, 0.f);

		for (int j = 0; j < m_NC[i]; j++) {
			std::unique_ptr<Fluid>& pcurr = fluidPs.at(m_Neighbor[i][j]);

			// Spiky Kernel
			glm::vec3 r = (p->predictPos - pcurr->predictPos) * m_Param[SPH_SIMSCALE];;
			r *= m_SpikyKern;			// grad of spiky (no minus)
			r *= (m_Param[SPH_SMOOTHRADIUS] - m_NDist[i][j]) * (m_Param[SPH_SMOOTHRADIUS] - m_NDist[i][j]) /
				m_NDist[i][j];
			r /= m_Param[SPH_RESTDENSITY];
			// End Spiky Kernel

			r *= (p->lambda + pcurr->lambda); // TODO scorr here
			// TODO scorr here

			p->deltaPos += r;
		}
	}
}

void FluidSystem::SPH_ApplyCorrections() {
	for (std::unique_ptr<Fluid>& p : fluidPs) {
		p->predictPos += p->deltaPos;

		// Perform collision detect and response
		if (p->predictPos.y < m_Vec[SPH_VOLMIN].y) { p->vel.y = 0.0; p->predictPos.y = m_Vec[SPH_VOLMIN].y + 0.01; }
		if (p->predictPos.y > m_Vec[SPH_VOLMAX].y) { p->vel.y = 0.0; p->predictPos.y = m_Vec[SPH_VOLMAX].y - 0.01; }

		if (p->predictPos.z < m_Vec[SPH_VOLMIN].z) { p->vel.z = 0.0; p->predictPos.z = m_Vec[SPH_VOLMIN].z + 0.01; }
		if (p->predictPos.z > m_Vec[SPH_VOLMAX].z) { p->vel.z = 0.0; p->predictPos.z = m_Vec[SPH_VOLMAX].z - 0.01; }

		if (p->predictPos.x < m_Vec[SPH_VOLMIN].x) { p->vel.x = 0.0; p->predictPos.x = m_Vec[SPH_VOLMIN].x + 0.01; }
		if (p->predictPos.x > m_Vec[SPH_VOLMAX].x) { p->vel.x = 0.0; p->predictPos.x = m_Vec[SPH_VOLMAX].x - 0.01; }
	}
}

void FluidSystem::Advance() {
	for (std::unique_ptr<Fluid>& p : fluidPs) {
		//std::cout << "p pos: " << p->pos.x << " " << p->pos.y << " " << p->pos.z << std::endl;
		//std::cout << "p pred: " << p->predictPos.x << " " << p->predictPos.y << " " << p->predictPos.z << std::endl;
		//std::cout << "p delP: " << p->deltaPos.x << " " << p->deltaPos.y << " " << p->deltaPos.z << std::endl;
		//std::cout << "p vel: " << p->vel.x << " " << p->vel.y << " " << p->vel.z << std::endl;
		//std::cout << "p dens: " << p->density << std::endl;
		//std::cout << "p lamb: " << p->lambda << std::endl << std::endl;

		p->vel = (p->predictPos - p->pos) / m_DT;
		//p->vel *= m_Param[SPH_SIMSCALE]; // ?? not sure if necesary
		//p->vel /= m_Param[SPH_SIMSCALE];
		// done
		//std::cout << glm::to_string(p->vel) << std::endl;
		//// apply vorticity confinement here?
		p->pos = p->predictPos;
	}
}

///--------------older methods

void FluidSystem::SPH_ComputeForceGridNC()
{
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		p->sph_force = glm::vec3(0, 0, 0);
		for (int j = 0; j < m_NC[i]; j++) {
			std::unique_ptr<Fluid>& pcurr = fluidPs.at(m_Neighbor[i][j]);
			float c = (m_Param[SPH_SMOOTHRADIUS] - m_NDist[i][j]);
			// m_NDist[i][j] == lenR, simscale distance
			float pterm = -0.5f * c * m_SpikyKern * (p->pressure + pcurr->pressure);
			float dterm = c * (1.f / p->density) * (1.f / pcurr->density);
			float vterm = m_LapKern * m_Param[SPH_VISC];
			p->sph_force += (pterm * (p->pos - pcurr->pos) +
				vterm * (pcurr->vel_eval - p->vel_eval)) * dterm;
		}
	}
}

void FluidSystem::AdvanceOld() {
	glm::vec3 norm, z;
	glm::vec3 dir, accel;
	glm::vec3 vnext;
	glm::vec3 min, max;
	double adj;
	float speed, diff;

	float stiff = m_Param[SPH_EXTSTIFF];
	float damp = m_Param[SPH_EXTDAMP];
	float radius = m_Param[SPH_PRADIUS];
	min = m_Vec[SPH_VOLMIN];
	max = m_Vec[SPH_VOLMAX];
	float ss = m_Param[SPH_SIMSCALE];

	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		// Compute Acceleration		
		accel = p->sph_force;
		accel *= m_Param[SPH_PMASS];

		// Velocity limiting 
		speed = accel.x * accel.x + accel.y * accel.y + accel.z * accel.z;
		if (speed > m_Param[SPH_LIMIT] * m_Param[SPH_LIMIT]) {
			accel *= m_Param[SPH_LIMIT] / sqrt(speed);
		}

		// Boundary Conditions
		// Z-axis walls (floor)
		diff = 2 * radius - (p->pos.z - min.z - (p->pos.x - m_Vec[SPH_VOLMIN].x) * m_Param[BOUND_ZMIN_SLOPE]) * ss;
		if (diff > EPSILON) {
			norm = glm::vec3(-m_Param[BOUND_ZMIN_SLOPE], 0, 1.0 - m_Param[BOUND_ZMIN_SLOPE]);
			adj = stiff * diff - damp * glm::dot(norm, p->vel_eval);
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}

		diff = 2 * radius - (max.z - p->pos.z) * ss;
		if (diff > EPSILON) {
			norm = glm::vec3(0, 0, -1);
			adj = stiff * diff - damp * glm::dot(norm, p->vel_eval);
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}

		// X-axis walls
		diff = 2 * radius - (p->pos.x - min.x + (sin(m_Time * 10.0) - 1 + (p->pos.y * 0.025) * 0.25) * m_Param[FORCE_XMIN_SIN]) * ss;
		//diff = 2 * radius - ( p->pos.x - min.x + (sin(m_Time*10.0)-1) * m_Param[FORCE_XMIN_SIN] )*ss;	
		if (diff > EPSILON) {
			norm = glm::vec3(1.0, 0, 0);
			adj = (m_Param[FORCE_XMIN_SIN] + 1) * stiff * diff - damp * glm::dot(norm, p->vel_eval);
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}

		diff = 2 * radius - (max.x - p->pos.x + (sin(m_Time * 10.0) - 1) * m_Param[FORCE_XMAX_SIN]) * ss;
		if (diff > EPSILON) {
			norm = glm::vec3(-1, 0, 0);
			adj = (m_Param[FORCE_XMAX_SIN] + 1) * stiff * diff - damp * glm::dot(norm, p->vel_eval);
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}

		// Y-axis walls
		diff = 2 * radius - (p->pos.y - min.y) * ss;
		if (diff > EPSILON) {
			norm = glm::vec3(0, 1, 0);
			adj = stiff * diff - damp * glm::dot(norm, p->vel_eval);
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}
		diff = 2 * radius - (max.y - p->pos.y) * ss;
		if (diff > EPSILON) {
			norm = glm::vec3(0, -1, 0);
			adj = stiff * diff - damp * glm::dot(norm, p->vel_eval);
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}

		// Plane gravity
		if (m_Param[PLANE_GRAV] > 0)
			accel += m_Vec[PLANE_GRAV_DIR];

		// Leapfrog Integration ----------------------------
		vnext = accel;
		vnext *= m_DT;
		vnext += p->vel;		// v(t+1/2) = v(t-1/2) + a(t) dt
		p->vel_eval = p->vel;
		p->vel_eval += vnext;
		p->vel_eval *= 0.5;		// v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5	  used to compute forces later
		p->vel = vnext;
		vnext *= m_DT / ss;
		p->pos += vnext;		// p(t+1) = p(t) + v(t+1/2) dt

		//if (p->pos.y < m_Vec[SPH_VOLMIN].y) { p->vel.y = 0.0; p->pos.y = m_Vec[SPH_VOLMIN].y + 0.01; }
		//if (p->pos.y > m_Vec[SPH_VOLMAX].y) { p->vel.y = 0.0; p->pos.y = m_Vec[SPH_VOLMAX].y - 0.01; }
		//if (p->pos.z < 0) { p->vel.z = 0.0; p->pos.z = 0.01; }
		//if (p->pos.z > m_Vec[SPH_VOLMAX].z) { p->vel.z = 0.0; p->pos.z = m_Vec[SPH_VOLMAX].z - 0.01; }
		//if (p->pos.x < m_Vec[SPH_VOLMIN].x) { p->vel.x = 0.0; p->pos.x = m_Vec[SPH_VOLMIN].x + 0.01; }
		//if (p->pos.x > m_Vec[SPH_VOLMAX].x) { p->vel.x = 0.0; p->pos.x = m_Vec[SPH_VOLMAX].x - 0.01; }
	}
	m_Time += m_DT;
}