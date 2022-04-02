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

#include "gl_helper.h"
#include <gl/glut.h>
#include <glm/gtx/norm.hpp>
#include <glm/gtx/string_cast.hpp>

#include "common_defs.h"
#include "fluid_system.h"

#define EPSILON			0.00001f			//for collision detection
#define GRAVITY glm::vec3(0, 0, -9.8)
#define SPH_RADIUS 0.1f
#define REST_DENSITY 6000.f

FluidSystem::FluidSystem()
{
}

void FluidSystem::Draw(float* view_mat, float rad)
{
	glEnable(GL_NORMALIZE);

	// point drawmode
	//glLoadMatrixf(view_mat);
	//dat = mBuf.data;
	//glBegin(GL_POINTS);
	//	for (auto& p : fluidPs) {
	//	glColor3f(RED(p->clr), GRN(p->clr), BLUE(p->clr));
	//	glVertex3f(p->pos.x, p->pos.y, p->pos.z);
	//	dat += mBuf.stride;
	//}
	//glEnd();
	glLoadMatrixf(view_mat);
	for (auto& f : fluidPs) {
		glm::vec3 scaledPos = f->pos;
		scaledPos /= 0.1;
		glPushMatrix();
		glTranslatef(scaledPos.x, scaledPos.y, scaledPos.z);
		glScalef(0.2, 0.2, 0.2);
		glColor4f(RED(f->clr), GRN(f->clr), BLUE(f->clr), ALPH(f->clr));
		drawSphere();
		glPopMatrix();
	}
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

void FluidSystem::SPH_CreateExample(int n, int nmax)
{
	m_Poly6Kern = 315.0f / (64.0f * 3.141592 * pow(SPH_RADIUS, 9));	// Wpoly6 kernel (denominator part) - 2003 Muller, p.4
	m_SpikyKern = -45.0f / (3.141592 * pow(SPH_RADIUS, 6));			// grad of spiky (no minus)

	m_Vec[SPH_VOLMIN] = glm::vec3(0, 0, 0);
	m_Vec[SPH_VOLMAX] = glm::vec3(10, 10, 10);
	// INIT MAX/ INIT MIN governs the two opposing corners of a square drop 
	//m_Vec[SPH_INITMIN] = glm::vec3(-5, -5, 10);
	//m_Vec[SPH_INITMAX] = glm::vec3(5, 5, 20);

	m_Vec[SPH_INITMIN] = m_Vec[SPH_VOLMIN];
	m_Vec[SPH_INITMAX] = m_Vec[SPH_VOLMAX]; //??????????

	gridSpaceDiag = (m_Vec[SPH_VOLMAX] - m_Vec[SPH_VOLMIN]); // / SPH_RADIUS;
	m_DT = 0.0083;

	for (int x = m_Vec[SPH_INITMIN].x; x < m_Vec[SPH_INITMAX].x; x++) {
		for (int y = m_Vec[SPH_INITMIN].y; y < m_Vec[SPH_INITMAX].y; y++) {
			for (int z = m_Vec[SPH_INITMIN].z; z < m_Vec[SPH_INITMAX].z; z++) {
				int ndx = fluidPs.size();
				fluidPs.push_back(std::make_unique<Fluid>());
				fluidPs.at(ndx)->pos = glm::vec3(x, y, z) * SPH_RADIUS * 0.5f; // * m_Param[SPH_RADIUS];
				fluidPs.at(ndx)->predictPos = glm::vec3(x, y, z);
				fluidPs.at(ndx)->clr = COLORA(0.5, 0.5, 0.5, 1);
			}
		}
	}

	// SETUP GRID
	totalGridCells = (int)(gridSpaceDiag.x * gridSpaceDiag.y * gridSpaceDiag.z);
	grid.reserve(totalGridCells);
	for (int _ = 0; _ < totalGridCells; ++_) {
		std::vector<int>& test = std::vector<int>();
		test.reserve(MAX_NEIGHBOR);
		grid.push_back(test);
	}

	for (int _ = 0; _ < fluidPs.size(); ++_) {
		neighbors.push_back(std::vector<int>());
	}
}

glm::vec3 FluidSystem::GetGridPos(glm::vec3 pos) {
	// modulo wrap not nedcessary?
	//glm::vec3 gridPos(int((pos.x - m_Param[SPH_VOLMIN].x) / SPH_RADIUS) % int(gridSpaceDiag.x),
	//				    int((pos.y - m_Param[SPH_VOLMIN].y) / SPH_RADIUS) % int(gridSpaceDiag.y),
	//				    int((pos.z - m_Param[SPH_VOLMIN].z) / SPH_RADIUS) % int(gridSpaceDiag.z));
	return glm::vec3(int(pos.x / SPH_RADIUS),
		int(pos.y / SPH_RADIUS),
		int(pos.z / SPH_RADIUS));
}

int FluidSystem::GetGridIndex(glm::vec3 gridPos) {
	// get index in grid space
	return gridPos.z * gridSpaceDiag.y * gridSpaceDiag.x +
		gridPos.y * gridSpaceDiag.x +
		gridPos.x;
}

void FluidSystem::PredictPositions() {
	for (std::unique_ptr<Fluid>& p : fluidPs) {
		glm::vec3 deltaVel = GRAVITY * m_DT; // a * dt = change in v

		// apply force to velocity (gravity)
		//p->vel += deltaVel;

		// consider limiting velocity
		//if (length vec > m_Param[SPH_LIMIT]) { p->vel = m_Param[SPH_LIMIT]; }
		p->predictPos = p->pos + (p->vel * m_DT);
	}
}

void FluidSystem::FindNeighbors() {
	for (int i = 0; i < totalGridCells; ++i) {
		grid.at(i).clear(); //gridcount is used to index into
	}

	//equivalinet of insertgrid / update grid finding the postns within the grid
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);

		glm::vec3 gridPos = GetGridPos(p->predictPos);
		int gIndex = GetGridIndex(gridPos);
		if (0 <= gIndex && gIndex < totalGridCells) {
			grid.at(gIndex).push_back(i); // maybe set a limit? (see MAX_NEIGHBOR)
		}
	}

	// equiv. Finding the Neighbors
	for (int i = 0; i < fluidPs.size(); ++i) {
		neighbors.at(i).clear(); // clear neighbors from prev. while were at it

		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		glm::vec3 gridPos = GetGridPos(p->predictPos);

		int SEARCH_SIZE = 1;
		// 2x2 neighborhood.
		for (int z = -SEARCH_SIZE; z <= SEARCH_SIZE; z++) {
		for (int y = -SEARCH_SIZE; y <= SEARCH_SIZE; y++) {
		for (int x = -SEARCH_SIZE; x <= SEARCH_SIZE; x++) {
			// go over each 2x2 of the GRID of the GRIDPOS of each fluid particle
			glm::vec3 n = gridPos + glm::vec3(x, y, z);
			if (n.x >= m_Vec[SPH_VOLMIN].x && n.x < m_Vec[SPH_VOLMAX].x &&
				n.y >= m_Vec[SPH_VOLMIN].y && n.y < m_Vec[SPH_VOLMAX].y &&
				n.z >= m_Vec[SPH_VOLMIN].z && n.z < m_Vec[SPH_VOLMAX].z) {
				int gIndex = GetGridIndex(n);
				for (int pIndex : grid.at(gIndex)) { // each 
					std::unique_ptr<Fluid>& pcurr = fluidPs.at(pIndex);
					float lenR = glm::length(p->predictPos - pcurr->predictPos);
					if (lenR <= SPH_RADIUS && lenR != 0) {
						neighbors.at(i).push_back(pIndex);
					}
				}
			}
		}}}
	}
}

void FluidSystem::Run()
{
	PredictPositions();
	FindNeighbors(); //?

	// PBF section
	for (int i = 0; i < 4; ++i) {
		ComputeDensity();
		ComputeLambda();
		ComputeCorrections();
		ApplyCorrections();
	}
	Advance();
}

void FluidSystem::ComputeDensity() {
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		p->density = 0.0;
		for (int j : neighbors.at(i)) { // for each neighbor
			std::unique_ptr<Fluid>& pcurr = fluidPs.at(j);

			// Poly Kernel
			float dist = glm::length(p->predictPos - pcurr->predictPos);
			float c = SPH_RADIUS * SPH_RADIUS - dist * dist;
			p->density += c * c * c;
		}
		p->density *= m_Poly6Kern;
	}
}

void FluidSystem::ComputeLambda() {
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		p->gradient = glm::vec3(0.0, 0.0, 0.0);
		float sumGradients = 0.0f;
		for (int j : neighbors.at(i)) { // for each neighbor
			std::unique_ptr<Fluid>& pcurr = fluidPs.at(j);

			// Spiky Kernel
			glm::vec3 r = (p->predictPos - pcurr->predictPos);
			float dist = glm::length(r);
			r *= m_SpikyKern;			// grad of spiky (no minus)
			r *= (SPH_RADIUS - dist) * (SPH_RADIUS - dist) / dist;
			r /= REST_DENSITY;
			// End Spiky Kernel

			sumGradients += glm::length2(r);
			p->gradient += r; // -= r; ?? - i think += b/c -45
		}
		sumGradients += glm::length2(p->gradient);
		float constraint = p->density / REST_DENSITY - 1.f; // real scale constraint
		p->lambda = -constraint / (sumGradients + 600.f); // maybe + 500 or so
	}
}

void FluidSystem::ComputeCorrections() {
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		p->deltaPos = glm::vec3(0.f, 0.f, 0.f);
		for (int j : neighbors.at(i)) { // for each neighbor
			std::unique_ptr<Fluid>& pcurr = fluidPs.at(j);
			glm::vec3 r = (p->predictPos - pcurr->predictPos);
			float dist = glm::length(r);
			//---------Calculate SCORR-----
			float cNum = SPH_RADIUS * SPH_RADIUS - dist * dist;
			float polyNum = cNum * cNum * cNum * m_Poly6Kern;

			// TODO, move magic numbers
			float delQ = 0.2f * SPH_RADIUS;
			float cDen = SPH_RADIUS * SPH_RADIUS - delQ * delQ;
			float polyDen = cDen * cDen * cDen * m_Poly6Kern;

			float frac = polyNum / polyDen;
			float K = 0.00001;
			float sCorr = -K * frac * frac * frac * frac;
			//------------End SCORR calculation-------

			// Spiky Kernel
			r *= (SPH_RADIUS - dist) * (SPH_RADIUS - dist) / dist;
			r *= m_SpikyKern;			// grad of spiky (no minus)
			// End Spiky Kernel

			//r *= (p->lambda + pcurr->lambda);
			r *= (p->lambda + pcurr->lambda + sCorr);

			p->deltaPos += r;
		}
		p->deltaPos /= REST_DENSITY; // (from kernel)
	}
}

void FluidSystem::ApplyCorrections() {
	for (std::unique_ptr<Fluid>& p : fluidPs) {
		p->predictPos += p->deltaPos;

		glm::vec3 test = m_Vec[SPH_VOLMIN];
		/*test.x += 10;*/
		glm::vec3 scaledMin = test * SPH_RADIUS;
		glm::vec3 scaledMax = m_Vec[SPH_VOLMAX] * SPH_RADIUS;
		// Perform collision detect and response
		if (p->predictPos.y < scaledMin.y) { p->vel.y = 0.0; p->predictPos.y = scaledMin.y + 0.01; }
		if (p->predictPos.y > scaledMax.y) { p->vel.y = 0.0; p->predictPos.y = scaledMax.y - 0.01; }

		if (p->predictPos.z < scaledMin.z) { p->vel.z = 0.0; p->predictPos.z = scaledMin.z + 0.01; }
		if (p->predictPos.z > scaledMax.z) { p->vel.z = 0.0; p->predictPos.z = scaledMax.z - 0.01; }

		if (p->predictPos.x < scaledMin.x) { p->vel.x = 0.0; p->predictPos.x = scaledMin.x + 0.01; }
		if (p->predictPos.x > scaledMax.x) { p->vel.x = 0.0; p->predictPos.x = scaledMax.x - 0.01; }
	}
}

void FluidSystem::Advance() {
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		p->vel = (p->predictPos - p->pos) / m_DT;
		//std::cout << "p pos: " << p->pos.x << " " << p->pos.y << " " << p->pos.z << std::endl;
		//std::cout << "p pred: " << p->predictPos.x << " " << p->predictPos.y << " " << p->predictPos.z << std::endl;
		//std::cout << "p delP: " << p->deltaPos.x << " " << p->deltaPos.y << " " << p->deltaPos.z << std::endl;
		//std::cout << "p vel: " << p->vel.x << " " << p->vel.y << " " << p->vel.z << std::endl;
		//std::cout << "p dens: " << p->density << std::endl;
		//std::cout << "p lamb: " << p->lambda << std::endl << std::endl;
	}

	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);

		// VISCOSITY
		glm::vec3 acc(0.f, 0.f, 0.f);
		for (int j : neighbors.at(i)) { // for each neighbor
			std::unique_ptr<Fluid>& pcurr = fluidPs.at(j);
			glm::vec3 vecDiff = pcurr->vel - p->vel;

			// check this, make sure its right!!
			float dist = glm::length(p->predictPos - pcurr->predictPos);
			float c = SPH_RADIUS * SPH_RADIUS - dist * dist;
			vecDiff *= c * c * c * m_Poly6Kern;

			acc += vecDiff;
		}
		float viscC = 0.01f;
		acc *= viscC;
		//std::cout << glm::to_string(p->vel) << " " << glm::to_string(acc) << std::endl;
		// Uncomment below in make it viscous
		p->vel += acc * m_DT;
		// END VISCOSITY

		p->pos = p->predictPos;
	}
}