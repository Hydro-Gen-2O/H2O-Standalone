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
	m_Param[MAX_FRAC] = 1.0;
	m_Param[POINT_GRAV] = 0.0;
	m_Param[PLANE_GRAV] = 1.0;

	m_Param[BOUND_ZMIN_SLOPE] = 0.0;
	m_Param[FORCE_XMAX_SIN] = 0.0;
	m_Param[FORCE_XMIN_SIN] = 0.0;
	m_Param[SPH_SMOOTHRADIUS] = 0.01;

	m_Vec[POINT_GRAV_POS].Set(0, 0, 50);
	m_Vec[PLANE_GRAV_DIR].Set(0, 0, -9.8);

	m_Param[SPH_SIMSCALE] = 0.004;			// unit size
	m_Param[SPH_VISC] = 0.2;			// pascal-second (Pa.s) = 1 kg m^-1 s^-1  (see wikipedia page on viscosity)
	m_Param[SPH_RESTDENSITY] = 600.0;			// kg / m^3
	m_Param[SPH_PMASS] = 0.00020543;		// kg
	m_Param[SPH_PRADIUS] = 0.004;			// m
	m_Param[SPH_INTSTIFF] = 1.00;
	m_Param[SPH_EXTSTIFF] = 10000.0;
	m_Param[SPH_EXTDAMP] = 256.0;
	m_Param[SPH_LIMIT] = 200.0;			// m / s

	// kernel computation
	m_Param[SPH_PDIST] = pow(m_Param[SPH_PMASS] / m_Param[SPH_RESTDENSITY], 1 / 3.0);
	m_Poly6Kern = 315.0f / (64.0f * 3.141592 * pow(m_Param[SPH_SMOOTHRADIUS], 9));	// Wpoly6 kernel (denominator part) - 2003 Muller, p.4
	m_SpikyKern = -45.0f / (3.141592 * pow(m_Param[SPH_SMOOTHRADIUS], 6));			// grad of spiky (no minus)
	m_LapKern = 45.0f / (3.141592 * pow(m_Param[SPH_SMOOTHRADIUS], 6));

	switch (n) {
	case -1:		// fluid drop
		m_Vec[SPH_VOLMIN].Set(-30, -30, 0);
		m_Vec[SPH_VOLMAX].Set(30, 30, 40);
		// INIT MAX/ INIT MIN governs the two opposing corners of a square drop 
		m_Vec[SPH_INITMIN].Set(-20, -26, 10);
		m_Vec[SPH_INITMAX].Set(20, 26, 40);
		break;
	case 0:		// Wave pool
		//-- TEST CASE: 2x2x2 grid, 32 particles.  NOTE: Set PRADIUS to 0.0004 to reduce wall influence
		m_Vec[SPH_VOLMIN].Set(-30, -30, 0);
		m_Vec[SPH_VOLMAX].Set(30, 30, 40);
		m_Vec[SPH_INITMIN].Set(-20, -26, 10);
		m_Vec[SPH_INITMAX].Set(20, 26, 40);

		m_Param[FORCE_XMIN_SIN] = 12.0;
		m_Param[BOUND_ZMIN_SLOPE] = 0.05;
		break;
	}

	m_DT = 0.003; //  0.001;			// .001 = for point grav
}

void FluidSystem::SPH_CreateExample(int n, int nmax)
{
	// setting up stuff
	SPH_Setup(n);
	fluidPs.clear();
	maxPoints = nmax;

	float ss = m_Param[SPH_PDIST] * 0.87 / m_Param[SPH_SIMSCALE];
	//printf ( "Spacing: %f\n", ss);
	Vector3DF min = m_Vec[SPH_INITMIN];
	Vector3DF max = m_Vec[SPH_INITMAX];
	float dx = max.x - min.x;
	float dy = max.y - min.y;
	float dz = max.z - min.z;
	for (float z = max.z; z >= min.z; z -= ss) {
		for (float y = min.y; y <= max.y; y += ss) {
			for (float x = min.x; x <= max.x; x += ss) {
				int ndx;
				if (fluidPs.size() < maxPoints - 2) {
					ndx = fluidPs.size();
					fluidPs.push_back(std::make_unique<Fluid>());
				}
				else {
					ndx = rand() % fluidPs.size();
				}
				fluidPs.at(ndx)->pos.Set(x, y, z);
				fluidPs.at(ndx)->clr = COLORA((x - min.x) / dx, (y - min.y) / dy, (z - min.z) / dz, 1);
			}
		}
	}
	Grid_Setup(
		m_Vec[SPH_VOLMIN],
		m_Vec[SPH_VOLMAX],
		m_Param[SPH_SIMSCALE],
		m_Param[SPH_SMOOTHRADIUS] * 2.0,
		1.0);
}

void FluidSystem::Run()
{
	Grid_InsertParticles();
	SPH_ComputeDensity(); // does the pressure comp for forces
	SPH_ComputeLambda();
	SPH_ComputeCorrections();

	//SPH_ComputeForceGridNC();
	Advance();
}

void FluidSystem::SPH_ComputeDensity()
{
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		float sum = 0.0;
		m_NC[i] = 0;
		p->gradient.Set(0.0, 0.0, 0.0);
		Grid_FindCells(p->pos, m_Param[SPH_SMOOTHRADIUS] / m_Param[SPH_SIMSCALE]);
		for (int cell = 0; cell < 8; cell++) {
			if (m_GridCell[cell] != -1) {
				int pndx = m_Grid[m_GridCell[cell]];
				while (pndx != -1) {
					std::unique_ptr<Fluid>& pcurr = fluidPs.at(pndx);
					if (pndx == i) { pndx = pcurr->next; continue; }
					float dx = (p->pos.x - pcurr->pos.x) * m_Param[SPH_SIMSCALE];	// dist in cm
					float dy = (p->pos.y - pcurr->pos.y) * m_Param[SPH_SIMSCALE];
					float dz = (p->pos.z - pcurr->pos.z) * m_Param[SPH_SIMSCALE];
					float dsq = (dx * dx + dy * dy + dz * dz);
					float lenR = sqrt(dsq);
					if (lenR < m_Param[SPH_SMOOTHRADIUS]) {
						float c = m_Param[SPH_SMOOTHRADIUS] * m_Param[SPH_SMOOTHRADIUS] - dsq;
						sum += c * c * c;

						//gradient calculation
						Vector3DF r(dx, dy, dz);
						r *= m_SpikyKern;
						r *= (m_Param[SPH_SMOOTHRADIUS] - lenR) * (m_Param[SPH_SMOOTHRADIUS] - lenR);
						std::cout << " LEN: " << lenR << std::endl;
						r /= lenR;
						r /= m_Param[SPH_RESTDENSITY];
						p->gradient += r;

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
		p->density = sum * m_Param[SPH_PMASS] * m_Poly6Kern;
		p->pressure = (p->density - m_Param[SPH_RESTDENSITY]) * m_Param[SPH_INTSTIFF];
	}
}

void FluidSystem::SPH_ComputeLambda() {
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		float constraint = p->density / m_Param[SPH_RESTDENSITY] - 1.f;
		float denominator = 0.f;

		for (int j = 0; j < m_NC[i]; j++) {
			std::unique_ptr<Fluid>& pcurr = fluidPs.at(m_Neighbor[i][j]);
			//denominator += (pcurr->gradient).Dot(pcurr->gradient);
			denominator += (pcurr->gradient).Length() * (pcurr->gradient).Length();
		}
		p->lambda = -constraint / (denominator + 500); // maybe + 500 or so
	}
}

void FluidSystem::SPH_ComputeCorrections() {
	//// the denominator of the surface correction (artificial pressure)
	//// and we choose k = 0.1, deltaq = 0.1h, n = 4
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		p->deltaPos.Set(0.f,0.f,0.f);

		for (int j = 0; j < m_NC[i]; j++) {
			std::unique_ptr<Fluid>& pcurr = fluidPs.at(m_Neighbor[i][j]);
			float magnitudeR = m_NDist[i][j];
			float dx = (p->pos.x - pcurr->pos.x) * m_Param[SPH_SIMSCALE];
			float dy = (p->pos.y - pcurr->pos.y) * m_Param[SPH_SIMSCALE];
			float dz = (p->pos.z - pcurr->pos.z) * m_Param[SPH_SIMSCALE];

			Vector3DF r(dx, dy, dz);
			r *= m_SpikyKern * // again, need negative?
				(m_Param[SPH_SMOOTHRADIUS] - magnitudeR) * (m_Param[SPH_SMOOTHRADIUS] - magnitudeR)
				/ magnitudeR *
				(p->lambda + pcurr->lambda) / 
				m_Param[SPH_RESTDENSITY];
			p->deltaPos += r;
		}
	}
}

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

void FluidSystem::Advance()
{
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		std::cout << "p delP: " << p->deltaPos.x << " " << p->deltaPos.y << " " << p->deltaPos.z << std::endl;
		std::cout << "p grad: " << p->gradient.x << " " << p->gradient.y << " " << p->gradient.z << std::endl;
		std::cout << "p dens: " << p->density << std::endl;
		std::cout << "p lamb: " << p->lambda << std::endl << std::endl;

		//Vector3DF tmp = p->pos;
		p->pos += p->deltaPos;

		//p->vel = p->pos;
		//p->vel -= tmp;
		//p->vel /= m_DT;

		//// apply vorticity confinement here?

		//p->pos += p->vel;
		//p->pos *= m_DT;
	}

	//Vector3DF norm, z;
	//Vector3DF dir, accel;
	//Vector3DF vnext;
	//Vector3DF min, max;
	//double adj;
	//float ss, radius;
	//float stiff, damp, speed, diff;

	//stiff = m_Param[SPH_EXTSTIFF];
	//damp = m_Param[SPH_EXTDAMP];
	//radius = m_Param[SPH_PRADIUS];
	//min = m_Vec[SPH_VOLMIN];
	//max = m_Vec[SPH_VOLMAX];
	//ss = m_Param[SPH_SIMSCALE];

	//for (int i = 0; i < fluidPs.size(); ++i) {
	//	std::unique_ptr<Fluid>& p = fluidPs.at(i);
	//	// Compute Acceleration		
	//	accel = p->sph_force;
	//	accel *= m_Param[SPH_PMASS];

	//	// Velocity limiting 
	//	speed = accel.x * accel.x + accel.y * accel.y + accel.z * accel.z;
	//	if (speed > m_Param[SPH_LIMIT] * m_Param[SPH_LIMIT]) {
	//		accel *= m_Param[SPH_LIMIT] / sqrt(speed);
	//	}

	//	// Boundary Conditions
	//	// Z-axis walls (floor)
	//	diff = 2 * radius - (p->pos.z - min.z - (p->pos.x - m_Vec[SPH_VOLMIN].x) * m_Param[BOUND_ZMIN_SLOPE]) * ss;
	//	if (diff > EPSILON) {
	//		norm.Set(-m_Param[BOUND_ZMIN_SLOPE], 0, 1.0 - m_Param[BOUND_ZMIN_SLOPE]);
	//		adj = stiff * diff - damp * norm.Dot(p->vel_eval);
	//		accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
	//	}

	//	diff = 2 * radius - (max.z - p->pos.z) * ss;
	//	if (diff > EPSILON) {
	//		norm.Set(0, 0, -1);
	//		adj = stiff * diff - damp * norm.Dot(p->vel_eval);
	//		accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
	//	}

	//	// X-axis walls
	//	diff = 2 * radius - (p->pos.x - min.x + (sin(m_Time * 10.0) - 1 + (p->pos.y * 0.025) * 0.25) * m_Param[FORCE_XMIN_SIN]) * ss;
	//	//diff = 2 * radius - ( p->pos.x - min.x + (sin(m_Time*10.0)-1) * m_Param[FORCE_XMIN_SIN] )*ss;	
	//	if (diff > EPSILON) {
	//		norm.Set(1.0, 0, 0);
	//		adj = (m_Param[FORCE_XMIN_SIN] + 1) * stiff * diff - damp * norm.Dot(p->vel_eval);
	//		accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
	//	}

	//	diff = 2 * radius - (max.x - p->pos.x + (sin(m_Time * 10.0) - 1) * m_Param[FORCE_XMAX_SIN]) * ss;
	//	if (diff > EPSILON) {
	//		norm.Set(-1, 0, 0);
	//		adj = (m_Param[FORCE_XMAX_SIN] + 1) * stiff * diff - damp * norm.Dot(p->vel_eval);
	//		accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
	//	}

	//	// Y-axis walls
	//	diff = 2 * radius - (p->pos.y - min.y) * ss;
	//	if (diff > EPSILON) {
	//		norm.Set(0, 1, 0);
	//		adj = stiff * diff - damp * norm.Dot(p->vel_eval);
	//		accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
	//	}
	//	diff = 2 * radius - (max.y - p->pos.y) * ss;
	//	if (diff > EPSILON) {
	//		norm.Set(0, -1, 0);
	//		adj = stiff * diff - damp * norm.Dot(p->vel_eval);
	//		accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
	//	}

	//	// Plane gravity
	//	if (m_Param[PLANE_GRAV] > 0)
	//		accel += m_Vec[PLANE_GRAV_DIR];

	//	// Point gravity
	//	if (m_Param[POINT_GRAV] > 0) {
	//		norm.x = (p->pos.x - m_Vec[POINT_GRAV_POS].x);
	//		norm.y = (p->pos.y - m_Vec[POINT_GRAV_POS].y);
	//		norm.z = (p->pos.z - m_Vec[POINT_GRAV_POS].z);
	//		norm.Normalize();
	//		norm *= m_Param[POINT_GRAV];
	//		accel -= norm;
	//	}

	//	// Leapfrog Integration ----------------------------
	//	vnext = accel;
	//	vnext *= m_DT;
	//	vnext += p->vel;		// v(t+1/2) = v(t-1/2) + a(t) dt
	//	p->vel_eval = p->vel;
	//	p->vel_eval += vnext;
	//	p->vel_eval *= 0.5;		// v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5	  used to compute forces later
	//	p->vel = vnext;
	//	vnext *= m_DT / ss;
	//	p->pos += vnext;		// p(t+1) = p(t) + v(t+1/2) dt

	//	if (m_Param[CLR_MODE] == 1.0) {
	//		adj = fabs(vnext.x) + fabs(vnext.y) + fabs(vnext.z) / 7000.0;
	//		adj = (adj > 1.0) ? 1.0 : adj;
	//		p->clr = COLORA(0, adj, adj, 1);
	//	}
	//	if (m_Param[CLR_MODE] == 2.0) {
	//		float v = 0.5 + (p->pressure / 1500.0);
	//		if (v < 0.1) v = 0.1;
	//		if (v > 1.0) v = 1.0;
	//		p->clr = COLORA(v, 1 - v, 0, 1);
	//	}
	//}
	//m_Time += m_DT;
}

//------------------------------------------------------ SPH Setup 
//
//  Range = +/- 10.0 * 0.006 (r) =	   0.12			m (= 120 mm = 4.7 inch)
//  Container Volume (Vc) =			   0.001728		m^3
//  Rest Density (D) =				1000.0			kg / m^3
//  Particle Mass (Pm) =			   0.00020543	kg						(mass = vol * density)
//  Number of Particles (N) =		4000.0
//  Water Mass (M) =				   0.821		kg (= 821 grams)
//  Water Volume (V) =				   0.000821     m^3 (= 3.4 cups, .21 gals)
//  Smoothing Radius (R) =             0.02			m (= 20 mm = ~3/4 inch)
//  Particle Radius (Pr) =			   0.00366		m (= 4 mm  = ~1/8 inch)
//  Particle Volume (Pv) =			   2.054e-7		m^3	(= .268 milliliters)
//  Rest Distance (Pd) =			   0.0059		m
//
//  Given: D, Pm, N
//    Pv = Pm / D			0.00020543 kg / 1000 kg/m^3 = 2.054e-7 m^3	
//    Pv = 4/3*pi*Pr^3    cuberoot( 2.054e-7 m^3 * 3/(4pi) ) = 0.00366 m
//     M = Pm * N			0.00020543 kg * 4000.0 = 0.821 kg		
//     V =  M / D              0.821 kg / 1000 kg/m^3 = 0.000821 m^3
//     V = Pv * N			 2.054e-7 m^3 * 4000 = 0.000821 m^3
//    Pd = cuberoot(Pm/D)    cuberoot(0.00020543/1000) = 0.0059 m 
//
// Ideal grid cell size (gs) = 2 * smoothing radius = 0.02*2 = 0.04
// Ideal domain size = k*gs/d = k*0.02*2/0.005 = k*8 = {8, 16, 24, 32, 40, 48, ..}
//    (k = number of cells, gs = cell size, d = simulation scale


//---------



// Compute Pressures - Using spatial grid, and also create neighbor table
void FluidSystem::SPH_ComputePressureGrid()
{
	
}

void FluidSystem::SPH_ComputeVorticityAndViscosity()
{
}