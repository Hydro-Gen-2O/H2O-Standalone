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

#include "common_defs.h"
#include "fluid_system.h"

#define EPSILON			0.00001f			//for collision detection

FluidSystem::FluidSystem()
{
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

void FluidSystem::SPH_Setup(int n) {
	//m_Param[BOUND_ZMIN_SLOPE] = 0.0;
	//m_Param[FORCE_XMAX_SIN] = 0.0;
	//m_Param[FORCE_XMIN_SIN] = 0.0;
	m_Param[SPH_INTSTIFF] = 0.50;
	m_Param[PLANE_GRAV] = 1.0;
	m_Param[SPH_VISC] = 0.2;			// pascal-second (Pa.s) = 1 kg m^-1 s^-1  (see wikipedia page on viscosity)

	m_Param[SPH_SMOOTHRADIUS] = 0.01; // 0.01;

	m_Vec[POINT_GRAV_POS].Set(0, 0, 50);
	m_Vec[PLANE_GRAV_DIR].Set(0, 0, -9.8);

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
		m_Vec[SPH_VOLMIN].Set(-10, -10, 0);
		m_Vec[SPH_VOLMAX].Set(10, 10, 30);
		// INIT MAX/ INIT MIN governs the two opposing corners of a square drop 
		m_Vec[SPH_INITMIN].Set(-5, -5, 10);
		m_Vec[SPH_INITMAX].Set(5, 5, 20);
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

Vector3DF FluidSystem::spikyKern(Vector3DF pos1, Vector3DF pos2) {
	float dx = (pos1.x - pos2.x) * m_Param[SPH_SIMSCALE];
	float dy = (pos1.y - pos2.y) * m_Param[SPH_SIMSCALE];
	float dz = (pos1.z - pos2.z) * m_Param[SPH_SIMSCALE];
	float dsq = (dx * dx + dy * dy + dz * dz);
	float lenR = sqrt(dsq);
	// assume neighborhood table arleady checked, else uncomment below
	//if (lenR > m_Param[SPH_SMOOTHRADIUS] || lenR == 0) {
	//	return Vector3DF(0.0, 0.0, 0.0);
	//}

	Vector3DF r(dx, dy, dz);
	r *= -45.0f / (3.141592 * pow(m_Param[SPH_SMOOTHRADIUS], 6));			// grad of spiky (no minus)

	//??????????????????
	//r *= (m_Param[SPH_SMOOTHRADIUS] / m_Param[SPH_SIMSCALE] - lenR) *
	//	(m_Param[SPH_SMOOTHRADIUS] / m_Param[SPH_SIMSCALE] - lenR);

	r *= (m_Param[SPH_SMOOTHRADIUS] - lenR) * (m_Param[SPH_SMOOTHRADIUS] - lenR);
	r /= lenR;
	return r;
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

	Vector3DF min = m_Vec[SPH_INITMIN];
	Vector3DF max = m_Vec[SPH_INITMAX];
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
				fluidPs.at(ndx)->pos.Set(x, y, z);
				fluidPs.at(ndx)->predictPos.Set(x, y, z);
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
		Vector3DF deltaVel = m_Vec[PLANE_GRAV_DIR];
		//deltaVel *= m_Param[SPH_PMASS]; // m * a
		deltaVel *= m_DT; // a * dt = change in v

		// TODO TOCHANGE NOTE apply force to velocity (gravity)
		//p->vel += deltaVel;

		// limit velocity ? 
		//if (p->vel.Length() > m_Param[SPH_LIMIT]) {
		//	p->vel = m_Param[SPH_LIMIT];
		//}

		Vector3DF velDist = p->vel;
		velDist *= m_DT; // / m_Param[SPH_SIMSCALE]; // ?????????
		p->predictPos = p->pos;
		p->predictPos += velDist;
	}
}

void FluidSystem::SPH_FindNeighbors(bool PBF) {
	float dx, dy, dz;

	Grid_InsertParticles();
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		float sum = 0.0;
		m_NC[i] = 0;
		if (PBF) {
			Grid_FindCells(p->predictPos, m_Param[SPH_SMOOTHRADIUS] / m_Param[SPH_SIMSCALE]);
		}
		else {
			Grid_FindCells(p->pos, m_Param[SPH_SMOOTHRADIUS] / m_Param[SPH_SIMSCALE]);
		}
		for (int cell = 0; cell < 8; cell++) {
			if (m_GridCell[cell] != -1) {
				int pndx = m_Grid[m_GridCell[cell]];
				while (pndx != -1) {
					std::unique_ptr<Fluid>& pcurr = fluidPs.at(pndx);
					if (pndx == i) { pndx = pcurr->next; continue; }
					if (PBF) {
						dx = (p->predictPos.x - pcurr->predictPos.x) * m_Param[SPH_SIMSCALE];	// dist in cm
						dy = (p->predictPos.y - pcurr->predictPos.y) * m_Param[SPH_SIMSCALE];
						dz = (p->predictPos.z - pcurr->predictPos.z) * m_Param[SPH_SIMSCALE];
					}
					else {
						dx = (p->pos.x - pcurr->pos.x) * m_Param[SPH_SIMSCALE];	// dist in cm
						dy = (p->pos.y - pcurr->pos.y) * m_Param[SPH_SIMSCALE];
						dz = (p->pos.z - pcurr->pos.z) * m_Param[SPH_SIMSCALE];
					}
					float dsq = (dx * dx + dy * dy + dz * dz);
					float lenR = sqrt(dsq);
					if (lenR != 0 && lenR < m_Param[SPH_SMOOTHRADIUS]) {
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

void FluidSystem::SPH_ComputeDensity()
{
	float dx, dy, dz;
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
		p->pressure = (p->density - m_Param[SPH_RESTDENSITY]) *m_Param[SPH_INTSTIFF];
	}
}

void FluidSystem::SPH_ComputeLambda() {
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);

		//WHICH?
		//float constraint = p->density / m_Param[SPH_RESTDENSITY] - 1.f;
		float constraint = p->density / m_Param[SPH_RESTDENSITY] / m_Param[SPH_SIMSCALE] - 1.f;

		p->gradient.Set(0.0, 0.0, 0.0);
		float sumGradients = 0.0f;

		for (int j = 0; j < m_NC[i]; j++) {
			std::unique_ptr<Fluid>& pcurr = fluidPs.at(m_Neighbor[i][j]);
			//gradient calculation
			Vector3DF r = spikyKern(p->predictPos, pcurr->predictPos);
			r /= m_Param[SPH_RESTDENSITY];

			//????
			//r *= m_Param[SPH_SIMSCALE]; // ? helps for some reason
			//??

			sumGradients += (r.Length() * r.Length());
			p->gradient += r; // -= r; ?? - i think += b/c -45
		}
		sumGradients += p->gradient.Length() * p->gradient.Length();
		p->lambda = -constraint / (sumGradients); // maybe + 500 or so
	}
}

void FluidSystem::SPH_ComputeCorrections() {
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		p->deltaPos.Set(0.f, 0.f, 0.f);

		for (int j = 0; j < m_NC[i]; j++) {
			std::unique_ptr<Fluid>& pcurr = fluidPs.at(m_Neighbor[i][j]);
			Vector3DF r = spikyKern(p->predictPos, pcurr->predictPos);
			r /= m_Param[SPH_RESTDENSITY];
			r *= (p->lambda + pcurr->lambda); // TODO scorr here

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

		if (p->predictPos.z < 0) { p->vel.z = 0.0; p->predictPos.z = 0.01; }
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

		p->vel = p->predictPos;
		p->vel -= p->pos;
		p->vel /= m_DT;
		//p->vel *= m_Param[SPH_SIMSCALE]; // ?? not sure if necesary

		// done

		//// apply vorticity confinement here?
		p->pos = p->predictPos;
	}
}

///--------------older methods

void FluidSystem::SPH_ComputeForceGridNC()
{
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		p->sph_force.Set(0, 0, 0);
		for (int j = 0; j < m_NC[i]; j++) {
			std::unique_ptr<Fluid>& pcurr = fluidPs.at(m_Neighbor[i][j]);
			float dx = (p->pos.x - pcurr->pos.x) * m_Param[SPH_SIMSCALE];		// dist in cm
			float dy = (p->pos.y - pcurr->pos.y) * m_Param[SPH_SIMSCALE];
			float dz = (p->pos.z - pcurr->pos.z) * m_Param[SPH_SIMSCALE];
			float c = (m_Param[SPH_SMOOTHRADIUS] - m_NDist[i][j]);
			float pterm = -0.5f * c * m_SpikyKern * (p->pressure + pcurr->pressure) / m_NDist[i][j];
			float dterm = c * (1.f / p->density) * (1.f / pcurr->density);
			float vterm = m_LapKern * m_Param[SPH_VISC];
			p->sph_force.x += (pterm * dx + vterm * (pcurr->vel_eval.x - p->vel_eval.x)) * dterm;
			p->sph_force.y += (pterm * dy + vterm * (pcurr->vel_eval.y - p->vel_eval.y)) * dterm;
			p->sph_force.z += (pterm * dz + vterm * (pcurr->vel_eval.z - p->vel_eval.z)) * dterm;
		}
	}
}

void FluidSystem::AdvanceOld() {
	Vector3DF norm, z;
	Vector3DF dir, accel;
	Vector3DF vnext;
	Vector3DF min, max;
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
		//diff = 2 * radius - (p->pos.z - min.z - (p->pos.x - m_Vec[SPH_VOLMIN].x) * m_Param[BOUND_ZMIN_SLOPE]) * ss;
		//if (diff > EPSILON) {
		//	norm.Set(-m_Param[BOUND_ZMIN_SLOPE], 0, 1.0 - m_Param[BOUND_ZMIN_SLOPE]);
		//	adj = stiff * diff - damp * norm.Dot(p->vel_eval);
		//	accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		//}

		//diff = 2 * radius - (max.z - p->pos.z) * ss;
		//if (diff > EPSILON) {
		//	norm.Set(0, 0, -1);
		//	adj = stiff * diff - damp * norm.Dot(p->vel_eval);
		//	accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		//}

		//// X-axis walls
		//diff = 2 * radius - (p->pos.x - min.x + (sin(m_Time * 10.0) - 1 + (p->pos.y * 0.025) * 0.25) * m_Param[FORCE_XMIN_SIN]) * ss;
		////diff = 2 * radius - ( p->pos.x - min.x + (sin(m_Time*10.0)-1) * m_Param[FORCE_XMIN_SIN] )*ss;	
		//if (diff > EPSILON) {
		//	norm.Set(1.0, 0, 0);
		//	adj = (m_Param[FORCE_XMIN_SIN] + 1) * stiff * diff - damp * norm.Dot(p->vel_eval);
		//	accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		//}

		//diff = 2 * radius - (max.x - p->pos.x + (sin(m_Time * 10.0) - 1) * m_Param[FORCE_XMAX_SIN]) * ss;
		//if (diff > EPSILON) {
		//	norm.Set(-1, 0, 0);
		//	adj = (m_Param[FORCE_XMAX_SIN] + 1) * stiff * diff - damp * norm.Dot(p->vel_eval);
		//	accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		//}

		//// Y-axis walls
		//diff = 2 * radius - (p->pos.y - min.y) * ss;
		//if (diff > EPSILON) {
		//	norm.Set(0, 1, 0);
		//	adj = stiff * diff - damp * norm.Dot(p->vel_eval);
		//	accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		//}
		//diff = 2 * radius - (max.y - p->pos.y) * ss;
		//if (diff > EPSILON) {
		//	norm.Set(0, -1, 0);
		//	adj = stiff * diff - damp * norm.Dot(p->vel_eval);
		//	accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		//}

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

		if (m_Param[CLR_MODE] == 1.0) {
			adj = fabs(vnext.x) + fabs(vnext.y) + fabs(vnext.z) / 7000.0;
			adj = (adj > 1.0) ? 1.0 : adj;
			p->clr = COLORA(0, adj, adj, 1);
		}
		if (m_Param[CLR_MODE] == 2.0) {
			float v = 0.5 + (p->pressure / 1500.0);
			if (v < 0.1) v = 0.1;
			if (v > 1.0) v = 1.0;
			p->clr = COLORA(v, 1 - v, 0, 1);
		}
	}
	m_Time += m_DT;
}