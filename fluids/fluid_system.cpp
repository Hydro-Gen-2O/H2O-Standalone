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
	m_R2 = m_Param[SPH_SMOOTHRADIUS] * m_Param[SPH_SMOOTHRADIUS];
	m_Poly6Kern = 315.0f / (64.0f * 3.141592 * pow(m_Param[SPH_SMOOTHRADIUS], 9));	// Wpoly6 kernel (denominator part) - 2003 Muller, p.4
	m_SpikyKern = -45.0f / (3.141592 * pow(m_Param[SPH_SMOOTHRADIUS], 6));			// grad of spiky (no minus)
	m_LapKern = 45.0f / (3.141592 * pow(m_Param[SPH_SMOOTHRADIUS], 6));

	switch (n) {
	case -1:		// fluid drop
		m_Vec[SPH_VOLMIN].Set(-30, -30, 0);
		m_Vec[SPH_VOLMAX].Set(30, 30, 40);
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
	m_Param[SPH_SIMSIZE] = m_Param[SPH_SIMSCALE] * (m_Vec[SPH_VOLMAX].z - m_Vec[SPH_VOLMIN].z);

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
				int ndx = AddPointReuse();
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

int FluidSystem::AddPointReuse()
{
	xref ndx;
	if (fluidPs.size() < maxPoints - 2) {
		ndx = fluidPs.size();

		fluidPs.push_back(std::make_unique<Fluid>());
		fluidPs.back()->sph_force.Set(0, 0, 0);
		fluidPs.back()->vel.Set(0, 0, 0);
		fluidPs.back()->vel_eval.Set(0, 0, 0);
		fluidPs.back()->next = 0x0;
		fluidPs.back()->pressure = 0;
		fluidPs.back()->density = 0;
	}
	else {
		ndx = rand() % fluidPs.size();
	}
	return ndx;
}

void FluidSystem::Run()
{
	Grid_InsertParticles();
	SPH_ComputePressureGridOld();
	SPH_ComputeForceGridNC();
	Advance();
}

void FluidSystem::SPH_ComputePressureGridOld()
{
	float radius = m_Param[SPH_SMOOTHRADIUS] / m_Param[SPH_SIMSCALE];
	float mR2 = m_Param[SPH_SMOOTHRADIUS] * m_Param[SPH_SMOOTHRADIUS];

	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& vecP = fluidPs.at(i);
		float sum = 0.0;
		m_NC[i] = 0;
		Grid_FindCells(vecP->pos, radius);
		for (int cell = 0; cell < 8; cell++) {
			if (m_GridCell[cell] != -1) {
				int pndx = m_Grid[m_GridCell[cell]];
				while (pndx != -1) {
					std::unique_ptr<Fluid>& pcurr = fluidPs.at(pndx);
					if (pndx == i) { pndx = pcurr->next; continue; }
					float dx = (vecP->pos.x - pcurr->pos.x) * m_Param[SPH_SIMSCALE];	// dist in cm
					float dy = (vecP->pos.y - pcurr->pos.y) * m_Param[SPH_SIMSCALE];
					float dz = (vecP->pos.z - pcurr->pos.z) * m_Param[SPH_SIMSCALE];
					float dsq = (dx * dx + dy * dy + dz * dz);
					if (mR2 > dsq) {
						float c = m_R2 - dsq;
						sum += c * c * c;
						if (m_NC[i] < MAX_NEIGHBOR) {
							m_Neighbor[i][m_NC[i]] = pndx;
							m_NDist[i][m_NC[i]] = sqrt(dsq);
							m_NC[i]++;
						}
					}
					pndx = pcurr->next;
				}
			}
			m_GridCell[cell] = -1;
		}
		vecP->density = sum * m_Param[SPH_PMASS] * m_Poly6Kern;
		vecP->pressure = (vecP->density - m_Param[SPH_RESTDENSITY]) * m_Param[SPH_INTSTIFF];
		vecP->density = 1.0f / vecP->density;
	}
}


void FluidSystem::Advance()
{
	Vector3DF norm, z;
	Vector3DF dir, accel;
	Vector3DF vnext;
	Vector3DF min, max;
	double adj;
	float SL, SL2, ss, radius;
	float stiff, damp, speed, diff;
	SL = m_Param[SPH_LIMIT];
	SL2 = SL * SL;

	stiff = m_Param[SPH_EXTSTIFF];
	damp = m_Param[SPH_EXTDAMP];
	radius = m_Param[SPH_PRADIUS];
	min = m_Vec[SPH_VOLMIN];
	max = m_Vec[SPH_VOLMAX];
	ss = m_Param[SPH_SIMSCALE];

	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& vecP = fluidPs.at(i);
		// Compute Acceleration		
		accel = vecP->sph_force;
		accel *= m_Param[SPH_PMASS];

		// Velocity limiting 
		speed = accel.x * accel.x + accel.y * accel.y + accel.z * accel.z;
		if (speed > SL2) {
			accel *= SL / sqrt(speed);
		}

		// Boundary Conditions
		// Z-axis walls (floor)
		diff = 2 * radius - (vecP->pos.z - min.z - (vecP->pos.x - m_Vec[SPH_VOLMIN].x) * m_Param[BOUND_ZMIN_SLOPE]) * ss;
		if (diff > EPSILON) {
			norm.Set(-m_Param[BOUND_ZMIN_SLOPE], 0, 1.0 - m_Param[BOUND_ZMIN_SLOPE]);
			adj = stiff * diff - damp * norm.Dot(vecP->vel_eval);
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}

		diff = 2 * radius - (max.z - vecP->pos.z) * ss;
		if (diff > EPSILON) {
			norm.Set(0, 0, -1);
			adj = stiff * diff - damp * norm.Dot(vecP->vel_eval);
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}

		// X-axis walls
		diff = 2 * radius - (vecP->pos.x - min.x + (sin(m_Time * 10.0) - 1 + (vecP->pos.y * 0.025) * 0.25) * m_Param[FORCE_XMIN_SIN]) * ss;
		//diff = 2 * radius - ( p->pos.x - min.x + (sin(m_Time*10.0)-1) * m_Param[FORCE_XMIN_SIN] )*ss;	
		if (diff > EPSILON) {
			norm.Set(1.0, 0, 0);
			adj = (m_Param[FORCE_XMIN_SIN] + 1) * stiff * diff - damp * norm.Dot(vecP->vel_eval);
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}

		diff = 2 * radius - (max.x - vecP->pos.x + (sin(m_Time * 10.0) - 1) * m_Param[FORCE_XMAX_SIN]) * ss;
		if (diff > EPSILON) {
			norm.Set(-1, 0, 0);
			adj = (m_Param[FORCE_XMAX_SIN] + 1) * stiff * diff - damp * norm.Dot(vecP->vel_eval);
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}

		// Y-axis walls
		diff = 2 * radius - (vecP->pos.y - min.y) * ss;
		if (diff > EPSILON) {
			norm.Set(0, 1, 0);
			adj = stiff * diff - damp * norm.Dot(vecP->vel_eval);
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}
		diff = 2 * radius - (max.y - vecP->pos.y) * ss;
		if (diff > EPSILON) {
			norm.Set(0, -1, 0);
			adj = stiff * diff - damp * norm.Dot(vecP->vel_eval);
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}

		// Plane gravity
		if (m_Param[PLANE_GRAV] > 0)
			accel += m_Vec[PLANE_GRAV_DIR];

		// Point gravity
		if (m_Param[POINT_GRAV] > 0) {
			norm.x = (vecP->pos.x - m_Vec[POINT_GRAV_POS].x);
			norm.y = (vecP->pos.y - m_Vec[POINT_GRAV_POS].y);
			norm.z = (vecP->pos.z - m_Vec[POINT_GRAV_POS].z);
			norm.Normalize();
			norm *= m_Param[POINT_GRAV];
			accel -= norm;
		}

		// Leapfrog Integration ----------------------------
		vnext = accel;
		vnext *= m_DT;
		vnext += vecP->vel;		// v(t+1/2) = v(t-1/2) + a(t) dt
		vecP->vel_eval = vecP->vel;
		vecP->vel_eval += vnext;
		vecP->vel_eval *= 0.5;		// v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5	  used to compute forces later
		vecP->vel = vnext;
		vnext *= m_DT / ss;
		vecP->pos += vnext;		// p(t+1) = p(t) + v(t+1/2) dt

		if (m_Param[CLR_MODE] == 1.0) {
			adj = fabs(vnext.x) + fabs(vnext.y) + fabs(vnext.z) / 7000.0;
			adj = (adj > 1.0) ? 1.0 : adj;
			vecP->clr = COLORA(0, adj, adj, 1);
		}
		if (m_Param[CLR_MODE] == 2.0) {
			float v = 0.5 + (vecP->pressure / 1500.0);
			if (v < 0.1) v = 0.1;
			if (v > 1.0) v = 1.0;
			vecP->clr = COLORA(v, 1 - v, 0, 1);
		}
	}

	m_Time += m_DT;
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
//    (k = number of cells, gs = cell size, d = simulation scale)
void FluidSystem::SPH_ComputeForceGridNC()
{
	int i = 0;
	for (i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& vecP = fluidPs.at(i);
		Vector3DF force;
		force.Set(0, 0, 0);
		for (int j = 0; j < m_NC[i]; j++) {
			std::unique_ptr<Fluid>& vecPcurr = fluidPs.at(m_Neighbor[i][j]);
			float dx = (vecP->pos.x - vecPcurr->pos.x) * m_Param[SPH_SIMSCALE];		// dist in cm
			float dy = (vecP->pos.y - vecPcurr->pos.y) * m_Param[SPH_SIMSCALE];
			float dz = (vecP->pos.z - vecPcurr->pos.z) * m_Param[SPH_SIMSCALE];
			float c = (m_Param[SPH_SMOOTHRADIUS] - m_NDist[i][j]);
			float pterm = -0.5f * c * m_SpikyKern * (vecP->pressure + vecPcurr->pressure) / m_NDist[i][j];
			float dterm = c * vecP->density * vecPcurr->density;
			float vterm = m_LapKern * m_Param[SPH_VISC];
			force.x += (pterm * dx + vterm * (vecPcurr->vel_eval.x - vecP->vel_eval.x)) * dterm;
			force.y += (pterm * dy + vterm * (vecPcurr->vel_eval.y - vecP->vel_eval.y)) * dterm;
			force.z += (pterm * dz + vterm * (vecPcurr->vel_eval.z - vecP->vel_eval.z)) * dterm;
		}
		vecP->sph_force = force;
	}
}


//---------



// Compute Pressures - Using spatial grid, and also create neighbor table
void FluidSystem::SPH_ComputePressureGrid()
{
}

void FluidSystem::SPH_ComputeVorticityAndViscosity()
{}